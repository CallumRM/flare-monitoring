define  proc_library  21160110228x
" fm000     Katherine12  Ke
" drudg version 2019Sep25 compiled under FS 10.00.00-alpha2-66-g0a0fe4de-dirty
"< DBBC_DDC             rack >< FlexBuff recorder 1>
enddef
define  exper_initi   21160110228x
proc_library
sched_initi
mk5=dts_id?
mk5=os_rev?
mk5_status
dbbc3=version
dbbcstop
dbbcinit
ifinit
dbbcstart
fbuff
enddef
define  setup01       21160110344
"Recorder may be wired to vsi0 or vsi2
mk5c_mode=vdif,0x00003333,,64.000
mk5c_mode
enddef
define  dbbcstop      21160110228
dbbc3=core3h=1,stop
dbbc3=core3h=2,stop
dbbc3=core3h=3,stop
dbbc3=core3h=4,stop
dbbc3=core3h=5,stop
dbbc3=core3h=6,stop
!+2s
enddef
define  dbbcinit      21160110232
dbbc3=core3h=3,inputselect vsi1
dbbc3=core3h=4,inputselect vsi1
dbbc3=core3h=5,inputselect vsi1
dbbc3=core3h=6,inputselect vsi1
!+1s
dbbc3=core3h=3,vsi_samplerate 128000000 2
dbbc3=core3h=4,vsi_samplerate 128000000 2
dbbc3=core3h=5,vsi_samplerate 128000000 2
dbbc3=core3h=6,vsi_samplerate 128000000 2
!+1s
dbbc3=core3h=3,splitmode off
dbbc3=core3h=4,splitmode off
dbbc3=core3h=5,splitmode off
dbbc3=core3h=6,splitmode off
!+1s
dbbc3=core3h=3,vsi_bitmask 0x00003333
dbbc3=core3h=4,vsi_bitmask 0x00003333
dbbc3=core3h=5,vsi_bitmask 0x00003333
dbbc3=core3h=6,vsi_bitmask 0x00003333
!+2s
dbbc3=core3h=3,reset
dbbc3=core3h=4,reset
dbbc3=core3h=5,reset
dbbc3=core3h=6,reset
!+2s
dbbc3=core3h=3,vdif_station Ke
dbbc3=core3h=4,vdif_station Ke
dbbc3=core3h=5,vdif_station Ke
dbbc3=core3h=6,vdif_station Ke
dbbc3=core3h=3,vdif_frame 2 4 8000 ct=off
dbbc3=core3h=4,vdif_frame 2 4 8000 ct=off
dbbc3=core3h=5,vdif_frame 2 4 8000 ct=off
dbbc3=core3h=6,vdif_frame 2 4 8000 ct=off
dbbc3=core3h=3,regupdate vdif_header 3 131072 0x03ff0000
dbbc3=core3h=4,regupdate vdif_header 3 196608 0x03ff0000
dbbc3=core3h=5,regupdate vdif_header 3 262144 0x03ff0000
dbbc3=core3h=6,regupdate vdif_header 3 327680 0x03ff0000
!+1s
dbbc3=core3h=3,destination 0 192.168.1.9:46227
dbbc3=core3h=4,destination 0 192.168.1.10:46227
dbbc3=core3h=5,destination 0 192.168.1.11:46227
dbbc3=core3h=6,destination 0 192.168.1.12:46227
!+2s
dbbc3=core3h=3,timesync
!+3s
dbbc3=core3h=4,timesync
!+3s
dbbc3=core3h=5,timesync
!+3s
dbbc3=core3h=6,timesync
!+3s
dbbc3=pps_sync
!+2s
dbbc3=dbbc17=165.00,c,32.00
dbbc3=dbbc18=652.00,c,32.00
dbbc3=dbbc19=1267.00,c,32.00
dbbc3=dbbc20=1814.00,c,32.00
dbbc3=dbbc25=165.00,d,32.00
dbbc3=dbbc26=652.00,d,32.00
dbbc3=dbbc27=1267.00,d,32.00
dbbc3=dbbc28=1814.00,d,32.00
dbbc3=dbbc33=2448.00,e,32.00
dbbc3=dbbc34=2663.00,e,32.00
dbbc3=dbbc35=2713.00,e,32.00
dbbc3=dbbc36=2814.00,e,32.00
dbbc3=dbbc41=2448.00,f,32.00
dbbc3=dbbc42=2663.00,f,32.00
dbbc3=dbbc43=2713.00,f,32.00
dbbc3=dbbc44=2814.00,f,32.00
!+2s
enddef
define  ifinit        21160110338x
dbbc3=dbbcifc=2,agc,2
dbbc3=dbbcifd=2,agc,2
dbbc3=dbbcife=2,agc,2
dbbc3=dbbciff=2,agc,2
enddef
define  dbbcstart     21160110340x
dbbc3=core3h=3,start vdif
dbbc3=core3h=4,start vdif
dbbc3=core3h=5,start vdif
dbbc3=core3h=6,start vdif
enddef
define  fbuff         21160110343x
mk5=mode=vdif_8000-512-4-2
mk5=mode?
mk5=net_port=46227
mk5=net_port?
mk5=net_protocol=udpsnor:64M:256M:1
mk5=mtu=9000
mk5=mtu?
"adding datastreams - suppose jive5ab_datastream is in use
mk5=datastream=reset
mk5=datastream=clear
mk5=datastream=add:c:192.168.1.17/Ke.*
mk5=datastream=add:d:192.168.1.18/Ke.*
mk5=datastream=add:e:192.168.1.19/Ke.*
mk5=datastream=add:f:192.168.1.20/Ke.*
mk5=datastream?
pcaloff
enddef
