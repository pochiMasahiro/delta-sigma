v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
N 80 30 80 70 {lab=0}
N 130 -0 240 0 {lab=out}
N 130 -0 160 0 {lab=out}
N 160 0 160 10 {lab=out}
N 160 10 160 110 {lab=out}
N 30 10 30 110 {lab=out}
N 30 10 50 10 {lab=out}
N 30 110 160 110 {lab=out}
N -100 -10 -100 10 {lab=in}
N -40 -10 50 -10 {lab=in}
N -100 -10 -40 -10 {lab=in}
C {opamp.sym} 80 0 0 0 {name=x1}
C {devices/gnd.sym} 80 70 0 0 {name=l1 lab=0}
C {devices/lab_pin.sym} 240 0 0 0 {name=p1 sig_type=std_logic lab=out
}
C {devices/vsource.sym} -100 40 0 0 {name=V1 value="sin(2.5 1 1k)" savecurrent=false}
C {devices/gnd.sym} -100 70 0 0 {name=l2 lab=0}
C {devices/code_shown.sym} -350 0 0 0 {name=s1 only_toplevel=false value="
.options savecurrent
.control
save all
tran 1n 1u
plot v(out) v(in)
.endc
"}
C {devices/lab_pin.sym} -100 -10 0 0 {name=p2 sig_type=std_logic lab=in}
