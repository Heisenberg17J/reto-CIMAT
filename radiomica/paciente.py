from dataclasses import dataclass
from pathlib import Path
import nibabel as nib
import pandas as pd

MODALIDADES = ('t1', 't1ce', 't2', 'flair', 'seg')

@dataclass
class Paciente:
    paciente_id: str
    t1: Path
    t1ce: Path
    t2: Path 
    flair: Path
    seg: Path

    @classmethod
    def desde_fila(cls, fila):
        return cls(fila['paciente_id'], *(Path(fila[m]) for m in MODALIDADES))

    def rutas(self):
        return {m: getattr(self, m) for m in MODALIDADES}
    
    def cargar(self, modalidad='t1'):
        return nib.load(self.rutas()[modalidad])

    def affine(self, modalidad='t1'):
        return self.cargar(modalidad).affine
    
    def dim_mm(self, modalidad='t1'):
        return self.cargar(modalidad).header.get_zooms()
    
def cargar_pacientes(manifest='datos/brats2018/manifest.csv'):
    
    return [Paciente.desde_fila(f) for _, f in pd.read_csv(manifest).iterrows()]

