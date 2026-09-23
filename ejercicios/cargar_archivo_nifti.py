import nibabel as nib
import matplotlib.pyplot as plt


partes_cerebro = ["sagital", "coronal", "axial"]
 
# Cargar la imagen NIfTI
imagen = nib.load('data/Brats18_2013_2_1_flair.nii.gz')

# Extraer los datos numéricos }
data = imagen.get_fdata()


# Dimension del NIfTI (x, y, z)
print(data.shape)

indice_corte = int(input(f'digite indice de corte: '))


corte_2d = [ data[indice_corte, : , :] , data[:, indice_corte, :],  data[:, :, indice_corte] ]

#creamos 3 espacios para poner las imagenes
fig, axes = plt.subplots(1, 3, figsize=(10,10))


#Dibujamos las 3 imagenes
plt.title(f'Cortes de cerebro en corte - {indice_corte}')
for i in range(3):
    # Mostrar la imagen con Matplotlib
    axes[i].imshow(corte_2d[i].T, cmap='gray', origin='lower')
    axes[i].set_title(f'Corte NIfTI {partes_cerebro[i]}')
    axes[i].axis('off')  # Ocultar los ejes

plt.tight_layout()  # Ajusta los márgenes automáticamente
plt.show()







