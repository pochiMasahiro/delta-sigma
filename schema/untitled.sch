v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
N 10 20 30 20 {lab=inn}
N 10 -20 30 -20 {lab=inp}
N 70 50 90 50 {lab=xxx}
N 70 30 70 50 {lab=xxx}
N 70 -50 70 -30 {lab=out}
N 70 -50 90 -50 {lab=out}
C {devices/vcvs.sym} 70 0 0 0 {name=E1 value=1e6}
C {devices/ipin.sym} 10 -20 0 0 {name=p1 lab=inp}
C {devices/ipin.sym} 10 20 0 0 {name=p2 lab=inn
}
C {devices/iopin.sym} 90 -50 0 0 {name=p3 lab=out
}
C {devices/iopin.sym} 90 50 0 0 {name=p4 lab=xxx}
