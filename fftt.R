## Definición funcion FFT
fftt = function(t,x) {
N=length(x)
T=t[2]-t[1]
F=1/T
f = seq(0,F/2,F/N)
y=fft(x)
y=y/(N/2)
y[1]=y[1]/2
y=y[1:(floor(N/2)+1)]
  return(output<-list(a=f,b=y))
}