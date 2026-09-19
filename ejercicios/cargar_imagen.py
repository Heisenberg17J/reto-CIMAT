import cv2

#cargar imagen
imagen = cv2.imread("/data/Brats18_2013_2_1_t1.nii.gz")

cv2.imshow("cerebro", imagen)

cv2.waitKey(0)
cv2.destroyAllWindows()