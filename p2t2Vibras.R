#Tarea 2 - Diemen Delgado C.

library(optimbase)

T = 0.05999985969 #periodo de giro [s]
wt = 2*pi/T  #frecuencia natural de giro [rad/s]
t = seq(0, T, T/250)  #tiempo en el que se evaluan las funciones

#Representación del momento f como series de fourier
a0 = 7500/8

f = a0/2*ones(length(t), 1)


for (n in 1:500){
  an = 500*sin(15*n*wt*T/16)/(pi*n)
  bn = -500*(cos(15*n*wt*T/16)-1)/(pi*n)
  f = f+an*cos(n*wt*t)+bn*sin(n*wt*t)
}


#grafico de f
plot(t, f, type="l", xlab="Tiempo [s]", ylab="f(t) [N m]", 
     main = "Par de torsión bajo el cual está sometida la fresa")
grid(NULL, NULL)

#parametros del sistema

#rango de diametros, se fue acotando al ir iterando el siguiente algoritmo
rd = seq(0.032, 0.036, by=0.000001)  
G = 8e10  #modulo de corte
L = 0.5  #largo del eje
I = 0.1  #inercia de la fresa
thetaMax = 0.0174533  #theta=1°


#Se itera la respuesta del sistema en función del diametro 
for (d in rd){
  kt=G*pi*d^4/(8*L)
  wn=sqrt(kt/I)
  #construir respuesta
  theta = (a0/(2*kt))*ones(length(t), 1)
  #500 terminos de la serie de fourier
  for (n in 1:500){
    an = 500*sin(15*n*wt*T/16)/(pi*n)
    bn = -500*(cos(15*n*wt*T/16)-1)/(pi*n)
    thetacn = (an/I)/(wn^2-(n*wt)^2)*cos(n*wt*t)
    thetasn = (bn/I)/(wn^2-(n*wt)^2)*sin(n*wt*t)
    theta = theta + thetacn + thetasn
  }
  #se comprueba que la amplitud de la respuesta sea cercana a thetaMax
  #diametro válido si el error es menor a 0.01%
  if (max(-theta)>max(theta)){
    if (100*abs(max(-theta)-thetaMax)/thetaMax < 0.01){ 
      cat("Posible valor de d:", d, "[m], ")
    }
  }
  else{
    if((100*abs(max(theta)-thetaMax)/thetaMax < 0.01)){
      cat("Posible valor de d:", d, "[m], ")
    }
  }
}

