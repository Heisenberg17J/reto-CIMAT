import cv2

#cargar imagen
imagen = cv2.imread("cerebro.jpg")

cv2.imshow("cerebro", imagen)

cv2.waitKey(0)
cv2.destroyAllWindows()