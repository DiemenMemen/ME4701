# Tarea 3 - Diemen Delgado C.

# leer librerias
library(R.matlab)
library(crone)

# Se leen los datos de datos,mat
datos = readMat("C:/Users/dieme/Downloads/datos.mat") # cambiar al directorio correcto
str(datos)
t = datos$t  # tiempo [s]
ay = datos$a  # aceleración en y [m/s^2]
vy = datos$v  # velocidad en y[m/s]
yh = datos$d  # desplazamiento en y[m]

# Otros datos
m = 835  # masa del transformador [kg]
h = 1.52  # altura del centro de masas [m]
Iv = 0.003847  # modulo de flexión [m^3}
E = 365.42e9  # modulo de Young [Pa]
Ta = 5.5e6  # tensión admisible [Pa]



# Parte a: calculo de k y c
delta = log(0.192/0.00105)  # decremento logaritmico
n = 7  # numero de ciclos entre las mediciones de amplitud
xi = delta/sqrt(4*pi^2*n^2 + delta^2)  # razón de amortiguamiento
wd = 2*pi*2.46  # frecuencia natural amortiguada [rad/s]
wn = wd/sqrt(1-xi^2)  # frecuencia natural [rad/s]
k = wn^2*m  # rigidez del transformador [N/m]
cc = 2*m*wn  # amortiguamiento crítico [kg/s]
c = cc*xi  # amortiguamiento del transformador [kg/s]
"k=202285.6 [N/m], c=3056.9 [kg/s]"  # resultados


# Parte b: respuesta del transformador frente a un terremoto
N = length(t)  # numero de pasos
dt = t[N]/N  # paso de tiempo [s]

# Condiciones iniciales
xh = matrix(0, N, 1)  # posición en x [m]
vx = matrix(0, N, 1)  # velocidad en x [m/s]
ax = matrix(0, N, 1)  # aceleración en x [m/s^2]
f = matrix(0, N, 1)  # excitación en la base [N]

# se aplica el método de diferencias centrales
# constantes
a0 = 1/dt^2
a1 = 1/(2*dt)
a2 = 2*a0
a3 = 1/a2

xh0 =xh[1] - dt*vx[1] + a3*ax[1]
hm = a0*m + a1*c
invm = 1/hm

# primer paso de tiempo
f[1] = c*vy[1] + k*yh[1] 
hf = f[1] - (k -a2*m)*xh[1] - (a0*m - a1*c)*xh0
xh[2] = invm*hf

# pasos de tiempo siguientes
for (i in 3:N){
  f[i] = c*vy[i] + k*yh[i] 
  hf = f[i-1] - (k -a2*m)*xh[i-1] - (a0*m - a1*c)*xh[i-2]
  xh[i] = invm*hf
}

# gráfico de la respuesta xh
dev.new()
plot(t, xh, type="l" ,xlab ="Tiempo [s]", ylab ="Respuesta [m]",
     main="Respuesta del transformador durante un sismo" )
grid(NULL, NULL)

# cálculo de la velocidad en x
for (i in 2:N){
  vx[i] = (xh[i+1] - xh[i-1])/(2*dt)
}
vx[N] = vx[N-1]  # Se asume que la velocidad no varía mucho al final



# Parte c: se comienza calculando la fuerza horizontal fh [N]
fh = matrix(0, N, 1)
for (i in 1:N){
  fh[i] = c*(vx[i] - vy[i]) + k*(xh[i] - yh[i])
}
Mh = h*fh  # momento horizontal [N*m]
Eh = Mh/Iv  # esfuerzo horizontal [Pa]

# Se debe verificar que el esfuerzo horizontal máximo sea menor a Ta
EhMax = max(abs(Eh))  # valor absoluto del esfuerzo horizontal máximo [Pa] 
cat("|EhMax|=", EhMax*10^(-6), "[MPa]")  # |EhMax|=5.146676 [MPa]
if (EhMax < Ta) {
  print("El esfuerzo horizontal durante el sismo no supera el esfuerzo admisible de 5.5 MPa.")
}

