import cv2
import nibabel as nib

#cargar imagen
data = nib.load('data/Brats18_2013_2_1/Brats18_2013_2_1_flair.nii.gz').get_fdata()

centro_x = data.shape[0] // 2
centro_y = data.shape[1] // 2
centro_z = data.shape[2] // 2

corte_axial = data[:, : , centro_z]
corte_sagital = data[centro_x, : , :]
corte_coronal = data[:, centro_y , :]

def procesar_corte(corte_crudo):
    corte_rotado = cv2.rotate(corte_crudo ,cv2.ROTATE_90_COUNTERCLOCKWISE)

    corte_normalizado = cv2.normalize(
        corte_rotado,
        None,
        alpha=0,
        beta=255,
        norm_type=cv2.NORM_MINMAX,
        dtype=cv2.CV_8U
    )

    return corte_normalizado

img_sagital = procesar_corte(corte_sagital)


cv2.imshow("corte sagital",img_sagital)

cv2.waitKey(0) 
cv2.destroyAllWindows()