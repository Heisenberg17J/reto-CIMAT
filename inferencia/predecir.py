"""
INFERENCIA DE PUNTA A PUNTA: segmentacion + pronostico para un paciente nuevo
Reto CIMAT / BraTS 2018  -- PROTOTIPO DE INVESTIGACION, NO ES UNA HERRAMIENTA CLINICA

    4 resonancias en formato BraTS (t1, t1ce, t2, flair) + edad
     -> nnU-Net (5 modelos de los folds, promediados) -> posprocesado de ET (D24)
     -> z-score (D2) -> regiones WT/TC/ET (D3) -> PyRadiomics con el mismo YAML (D7-D9)
     -> Cox Elastic Net final (D30) -> riesgo, clase de supervivencia y dias estimados

Reutiliza las funciones del pipeline (radiomica/, segmentacion/, pronostico/), asi
que el paciente se procesa exactamente igual que los 285 del estudio.

Requisitos:
    - Entorno con PyRadiomics + PyTorch + nnU-Net (ver inferencia/README.md).
    - Pesos de nnU-Net en datos/nnunet_results/ (descargados de Drive), salvo que
      se pase --segmentacion con una mascara ya hecha (etiquetas BraTS 0/1/2/4).
    - modelos/pronostico_coxnet.joblib (pronostico/entrenar_final.py).

Entrada: una carpeta con los 4 NIfTI, llamados t1/t1ce/t2/flair.nii.gz o con el
nombre de BraTS (<id>_t1.nii.gz, ...). Formato BraTS: sin craneo, co-registradas,
voxel de 1 mm, 240 x 240 x 155.

Salida en --salida:
    segmentacion.nii.gz   etiquetas BraTS (1 necrosis, 2 edema, 4 realce)
    caracteristicas.csv   las 1146 caracteristicas radiomicas
    resultado.json        volumenes, riesgo, clase, dias y advertencias
    vista.png             cortes con la segmentacion superpuesta

Uso (desde la raiz del repo):
    python inferencia/predecir.py --caso CARPETA --edad 62
    python inferencia/predecir.py --caso CARPETA --edad 62 --segmentacion mascara.nii.gz
    python inferencia/predecir.py --caso CARPETA --segmentacion mascara.nii.gz   # sin edad: sin pronostico
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, str(max(1, (os.cpu_count() or 2) - 1)))

import argparse
import json
import logging
import shutil
import sys
import tempfile
from pathlib import Path

import joblib
import nibabel as nib
import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent.parent
for _d in ("radiomica", "segmentacion", "pronostico"):
    sys.path.insert(0, str(RAIZ / _d))

import brazos as b                                  # noqa: E402
import extraccion as e                              # noqa: E402
import normalizar as n                              # noqa: E402
import regiones as r                                # noqa: E402
from convertir_nnunet import NNUNET_A_BRATS, reasignar  # noqa: E402

MODALIDADES = ["t1", "t1ce", "t2", "flair"]
MODELO_NNUNET = Path("datos/nnunet_results/Dataset501_BraTS2018/"
                     "nnUNetTrainer_100epochs_ckpt5__nnUNetPlans__3d_fullres")
MODELO_PRONOSTICO = Path("modelos/pronostico_coxnet.joblib")
POSPROCESO = Path("resultados/segmentacion/postproceso_et.json")
CLASES = {0: "corta (< 300 dias, < 10 meses)", 1: "media (300-450 dias, 10-15 meses)",
          2: "larga (> 450 dias, > 15 meses)"}
AVISO = ("PROTOTIPO DE INVESTIGACION. No es una herramienta clinica ni sustituye el criterio "
         "medico. El pronostico depende casi solo de la edad (D27, D30).")


# ---------------------------------------------------------------------------
# 1. Entrada
# ---------------------------------------------------------------------------

def buscar_modalidades(carpeta):
    """Acepta t1.nii.gz o <id>_t1.nii.gz (ojo: t1 no debe confundirse con t1ce).
    Busca tambien un nivel mas abajo, como en segmentar_colab.ipynb."""
    if not carpeta.is_dir():
        raise SystemExit(f"--caso {carpeta} no existe o no es una carpeta. Debe ser la carpeta con las 4 "
                         f"resonancias del paciente (t1, t1ce, t2 y flair).")
    archivos = list(carpeta.glob("*.nii*")) + list(carpeta.glob("*/*.nii*"))
    rutas = {}
    for mod in MODALIDADES:
        candidatos = [f for f in archivos
                      if f.name.split(".")[0] == mod or f.name.split(".")[0].endswith(f"_{mod}")]
        if len(candidatos) != 1:
            vistos = [str(f.relative_to(carpeta)) for f in archivos] or "ninguno"
            raise SystemExit(f"en {carpeta} hay {len(candidatos)} archivos para '{mod}' (se espera 1). "
                             f"Archivos .nii encontrados: {vistos}. Deben llamarse {mod}.nii.gz o <id>_{mod}.nii.gz")
        rutas[mod] = candidatos[0]
    return rutas


def validar_formato(rutas):
    """Avisos si las imagenes no parecen estar en formato BraTS (no detiene el proceso)."""
    imgs = {m: nib.load(p) for m, p in rutas.items()}
    ref, avisos = imgs["t1"], []
    for m, im in imgs.items():
        if im.shape != ref.shape or not np.allclose(im.affine, ref.affine, atol=1e-3):
            avisos.append(f"{m}: forma o geometria distinta de t1 (no co-registradas)")
    if ref.shape != (240, 240, 155):
        avisos.append(f"tamano {ref.shape}; BraTS usa (240, 240, 155)")
    if not np.allclose(ref.header.get_zooms()[:3], 1, atol=1e-2):
        avisos.append(f"voxel {tuple(round(float(z), 2) for z in ref.header.get_zooms()[:3])} mm; BraTS usa 1 mm")
    fondo = float((np.asanyarray(ref.dataobj) == 0).mean())
    if fondo < 0.5:
        avisos.append(f"solo {fondo:.0%} de fondo en 0: probablemente conserva el craneo")
    return avisos


# ---------------------------------------------------------------------------
# 2. Segmentacion (nnU-Net + posprocesado de ET)
# ---------------------------------------------------------------------------

def registrar_entrenador():
    """nnU-Net necesita encontrar la clase del entrenador con la que se guardaron los pesos."""
    import nnunetv2
    destino = Path(nnunetv2.__path__[0]) / "training" / "nnUNetTrainer" / "variants" / \
        "nnUNetTrainer_100epochs_ckpt5.py"
    if not destino.exists():
        destino.write_text(
            "from nnunetv2.training.nnUNetTrainer.variants.training_length.nnUNetTrainer_Xepochs "
            "import nnUNetTrainer_100epochs\n\n\n"
            "class nnUNetTrainer_100epochs_ckpt5(nnUNetTrainer_100epochs):\n"
            "    def on_train_start(self):\n"
            "        super().on_train_start()\n"
            "        self.save_every = 5\n")


def segmentar_nnunet(rutas, trabajo, folds, tta, log):
    import torch
    from nnunetv2.inference.predict_from_raw_data import nnUNetPredictor

    if not (MODELO_NNUNET / "plans.json").exists():
        raise SystemExit(f"no encuentro los pesos de nnU-Net en {MODELO_NNUNET} (ver inferencia/README.md)")
    registrar_entrenador()
    disp = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log.info(f"nnU-Net en {disp} | folds {folds} | TTA {'si' if tta else 'no'}"
             + ("  (en CPU puede tardar varios minutos)" if disp.type == "cpu" else ""))
    pred = nnUNetPredictor(tile_step_size=0.5, use_gaussian=True, use_mirroring=tta,
                           device=disp, verbose=False, allow_tqdm=True)
    pred.initialize_from_trained_model_folder(str(MODELO_NNUNET), use_folds=tuple(folds),
                                              checkpoint_name="checkpoint_final.pth")
    pred.predict_from_files([[str(rutas[m]) for m in MODALIDADES]], [str(trabajo / "nnunet")],
                            save_probabilities=False, overwrite=True,
                            num_processes_preprocessing=1, num_processes_segmentation_export=1)
    img = nib.load(trabajo / "nnunet.nii.gz")
    return np.asanyarray(img.dataobj).astype(np.int16), img


def posprocesar_et(seg_nnunet, log):
    """D24: ET predicho por debajo del umbral -> necrosis/no realzado. Devuelve etiquetas BraTS."""
    umbral = json.loads(POSPROCESO.read_text())["umbral_casos_nuevos"]
    n_et = int((seg_nnunet == 3).sum())
    if 0 < n_et < umbral:
        seg_nnunet = seg_nnunet.copy()
        seg_nnunet[seg_nnunet == 3] = 2
        log.info(f"posprocesado: ET de {n_et} voxeles < {umbral} -> descartado")
    return reasignar(seg_nnunet, NNUNET_A_BRATS)


# ---------------------------------------------------------------------------
# 3. Radiomica (mismo pipeline que los bloques 3-5)
# ---------------------------------------------------------------------------

def extraer_caracteristicas(rutas, seg_brats, img_ref, trabajo, log):
    norm = {}
    for mod in MODALIDADES:
        img = nib.load(rutas[mod])
        z, _ = n.zscore(img.get_fdata())
        nueva = nib.Nifti1Image(z.astype(np.float32), img.affine, img.header)
        nueva.set_data_dtype(np.float32)
        norm[mod] = trabajo / f"norm_{mod}.nii.gz"
        nib.save(nueva, norm[mod])

    umbrales = r.leer_umbrales()
    tareas, estados, volumenes = [], {}, {}
    for region, etiquetas in r.REGIONES.items():
        mascara = np.isin(seg_brats, etiquetas)
        nvox, _, estado = r.estado_region(mascara, *umbrales)
        estados[region], volumenes[region] = estado, nvox / 1000
        if estado != "ok":
            log.info(f"region {region}: {estado} ({nvox} voxeles) -> sus caracteristicas quedan vacias")
            continue
        ruta = trabajo / f"mask_{region}.nii.gz"
        m_img = nib.Nifti1Image(mascara.astype(np.uint8), img_ref.affine, img_ref.header)
        m_img.set_data_dtype(np.uint8)
        nib.save(m_img, ruta)
        base = {"paciente_id": "caso", "region": region, "mascara": str(ruta)}
        tareas.append({**base, "tipo": "forma", "modalidad": None, "imagen": str(norm[r.MODALIDAD_FORMA])})
        tareas += [{**base, "tipo": "intensidad", "modalidad": mod, "imagen": str(norm[mod])} for mod in MODALIDADES]

    forma, intensidad = r.crear_extractores()
    caract, _ = e.extraer_paciente(pd.DataFrame(tareas), {"forma": forma, "intensidad": intensidad})
    return caract, estados, volumenes


# ---------------------------------------------------------------------------
# 4. Pronostico
# ---------------------------------------------------------------------------

def clase_desde_riesgo(modelo, riesgo):
    c = modelo["cortes"]
    return 0 if riesgo >= c["corta_si_riesgo_>="] else (1 if riesgo >= c["media_si_riesgo_>="] else 2)


def pronosticar(caract, edad, estados):
    modelo = joblib.load(MODELO_PRONOSTICO)
    faltan = [c for c in modelo["columnas_radiomicas"] if not np.isfinite(caract.get(c, np.nan))]
    if faltan:
        sin = sorted({c.split("_")[1] for c in faltan})
        return {"aplica": False,
                "motivo": f"faltan regiones {sin} (estado: {estados}). El modelo se entreno con HGG, todos "
                          f"con realce (ET), y no se puede aplicar sin el. Un ET vacio puede indicar un tumor "
                          f"de bajo grado, o un realce pequeno que el posprocesado descarto (D24)."}
    X = pd.DataFrame([{c: caract[c] for c in modelo["columnas_radiomicas"]}])
    if edad is None:
        # Sin edad no se predice: la edad lleva casi toda la senal (D30) e imputarla seria dar el promedio.
        # Se muestra, solo como referencia, como cambiaria el resultado segun la edad.
        ref = {}
        for e_ in (40, 50, 60, 70, 80):
            r_, d_ = b.predecir_coxnet(modelo, pd.DataFrame([{"edad": float(e_)}]), X)
            ref[str(e_)] = f"{CLASES[clase_desde_riesgo(modelo, r_[0])]}, ~{d_[0]:.0f} dias"
        return {"aplica": False,
                "motivo": "falta la edad. El modelo depende casi por completo de ella (D30) y no se imputa: "
                          "con la edad media, el resultado seria solo el promedio de la cohorte.",
                "solo_referencia_segun_edad (NO es una prediccion)": ref}
    clin = pd.DataFrame([{"edad": float(edad)}])
    riesgo, dias = b.predecir_coxnet(modelo, clin, X)
    clase = clase_desde_riesgo(modelo, riesgo[0])
    return {"aplica": True, "riesgo": float(riesgo[0]), "clase": CLASES[clase],
            "dias_mediana_estimada": float(dias[0]),
            "desempeno_esperado": modelo["meta"]["desempeno_esperado_cv_5x10 (D27)"]}


# ---------------------------------------------------------------------------
# 5. Figura
# ---------------------------------------------------------------------------

def vista(rutas, seg, destino, titulo):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    colores = ListedColormap([(0, 0, 0, 0), (0.9, 0.1, 0.1, 1), (0.1, 0.8, 0.2, 1), (1.0, 0.85, 0.0, 1)])
    s = np.select([seg == 1, seg == 2, seg == 4], [1, 2, 3], 0)     # 1 necrosis, 2 edema, 3 realce
    flair = np.asanyarray(nib.load(rutas["flair"]).dataobj)
    # Los tres planos pasan por el centro del nucleo (TC); si no hay, por el del tumor completo
    nucleo = np.isin(s, (1, 3)) if np.isin(s, (1, 3)).any() else s > 0
    cx, cy, cz = (np.round(np.argwhere(nucleo).mean(axis=0)).astype(int) if nucleo.any()
                  else np.array(s.shape) // 2)
    cortes = [("axial", lambda a: np.rot90(a[:, :, cz])), ("coronal", lambda a: np.rot90(a[:, cy, :])),
              ("sagital", lambda a: np.rot90(a[cx, :, :]))]
    fig, ejes = plt.subplots(1, 3, figsize=(13, 4.8))
    for eje, (nombre, corte) in zip(ejes, cortes):
        eje.imshow(corte(flair), cmap="gray")
        eje.imshow(corte(s), cmap=colores, vmin=0, vmax=3, alpha=0.45, interpolation="nearest")
        eje.set_title(f"FLAIR, {nombre}")
        eje.axis("off")
    fig.legend(handles=[Patch(color=colores(1), label="necrosis / no realzado"),
                        Patch(color=colores(2), label="edema"), Patch(color=colores(3), label="realce (ET)")],
               loc="lower center", ncol=3, frameon=False)
    fig.suptitle(titulo, fontsize=10)
    fig.tight_layout(rect=(0, 0.07, 1, 0.94))
    fig.savefig(destino, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--caso", required=True, type=Path, help="carpeta con t1, t1ce, t2 y flair")
    ap.add_argument("--edad", type=float, help="edad en anos; sin ella no hay pronostico (el modelo depende de ella)")
    ap.add_argument("--segmentacion", type=Path, help="mascara ya hecha (BraTS 0/1/2/4): omite nnU-Net")
    ap.add_argument("--folds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--tta", action="store_true", help="aumento en prueba (espejos): mejor y ~8x mas lento")
    ap.add_argument("--salida", type=Path, default=None)
    args = ap.parse_args()

    caso = args.caso.resolve()
    seg_ext = args.segmentacion.resolve() if args.segmentacion else None
    salida = (args.salida or Path("resultados_inferencia") / caso.name).resolve()
    os.chdir(RAIZ)

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    log = logging.getLogger("inferencia")
    import radiomics
    radiomics.setVerbosity(logging.ERROR)
    for ruidoso in ("radiomics", "pykwalify"):      # no propagar sus mensajes al log principal
        logging.getLogger(ruidoso).setLevel(logging.ERROR)
    log.info(AVISO + "\n")

    rutas = buscar_modalidades(caso)
    if seg_ext and not seg_ext.is_file():
        raise SystemExit(f"--segmentacion {seg_ext} no existe")
    salida.mkdir(parents=True, exist_ok=True)
    avisos = validar_formato(rutas)
    for a in avisos:
        log.warning(f"AVISO de formato: {a}")

    trabajo = Path(tempfile.mkdtemp(prefix="inferencia_"))
    try:
        if seg_ext:
            img_ref = nib.load(seg_ext)
            seg = np.asanyarray(img_ref.dataobj).astype(np.int16)
            log.info(f"segmentacion: {seg_ext} (sin nnU-Net)")
            fuente = f"externa: {seg_ext.name}"
        else:
            seg_nn, img_ref = segmentar_nnunet(rutas, trabajo, args.folds, args.tta, log)
            seg = posprocesar_et(seg_nn, log)
            fuente = f"nnU-Net, folds {args.folds}, TTA {'si' if args.tta else 'no'}"
        nueva = nib.Nifti1Image(seg.astype(np.uint8), img_ref.affine, img_ref.header)
        nueva.set_data_dtype(np.uint8)
        nib.save(nueva, salida / "segmentacion.nii.gz")

        log.info("radiomica: normalizacion, regiones y PyRadiomics...")
        caract, estados, vols = extraer_caracteristicas(rutas, seg, img_ref, trabajo, log)
        pd.DataFrame([caract]).to_csv(salida / "caracteristicas.csv", index=False)
    finally:
        shutil.rmtree(trabajo, ignore_errors=True)

    pron = pronosticar(caract, args.edad, estados)
    resultado = {"aviso": AVISO, "caso": str(caso), "edad": args.edad, "segmentacion": fuente,
                 "volumen_cm3": {k: round(v, 2) for k, v in vols.items()}, "estado_regiones": estados,
                 "avisos_formato": avisos, "pronostico": pron}
    (salida / "resultado.json").write_text(json.dumps(resultado, indent=2, ensure_ascii=False))

    titulo = (f"WT {vols['WT']:.1f} cm³ · TC {vols['TC']:.1f} cm³ · ET {vols['ET']:.1f} cm³ | "
              + (f"supervivencia {pron['clase']}, ~{pron['dias_mediana_estimada']:.0f} días"
                 if pron["aplica"] else "pronóstico: no aplica") + " | prototipo de investigación")
    vista(rutas, seg, salida / "vista.png", titulo)

    log.info(f"\nvolumen (cm³): {resultado['volumen_cm3']}")
    if pron["aplica"]:
        log.info(f"pronostico: {pron['clase']} | mediana estimada ~{pron['dias_mediana_estimada']:.0f} dias "
                 f"| riesgo {pron['riesgo']:+.3f}")
        log.info(f"desempeno esperado del modelo (CV): {pron['desempeno_esperado']}")
    else:
        log.info(f"pronostico: no aplica. {pron['motivo']}")
        for e_, txt in pron.get("solo_referencia_segun_edad (NO es una prediccion)", {}).items():
            log.info(f"   referencia, si tuviera {e_} anos: {txt}")
    log.info(f"\nsalida: {salida}/  (segmentacion.nii.gz, caracteristicas.csv, resultado.json, vista.png)")


if __name__ == "__main__":
    main()
