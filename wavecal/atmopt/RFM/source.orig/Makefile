
# Name of executable and main source
EXENAME  = rfm

OPTIM = -O3

#FF = ifort -132 *.for
FF = gfortran *.for

FLAGS   = $(OPTIM) 

COMMONS = adjcom.inc atmcom.inc chidat.inc co2dat.inc crvcom.inc cyccom.inc fincom.inc flgcom.inc fovcom.inc gascom.inc gflcom.inc gracom.inc grdcom.inc h2ockd.inc h2omtc.inc hflcom.inc hitcom.inc idxcon.inc ilscom.inc jaccom.inc jdxcon.inc lflcom.inc logcom.inc lutcom.inc mixdat.inc ntecom.inc obscom.inc optdat.inc outcom.inc phycon.inc pthcom.inc qalcom.inc qfncom.inc reclen.inc rejcom.inc rfmcon.inc rfmfil.inc rfmsiz.inc sfccom.inc shpcon.inc spccom.inc tabcom.inc tancom.inc tpscom.inc widcom.inc xflcom.inc xsccom.inc


rfm: $(COMMONS)  Makefile
	$(FF) $(LFLAGS) -o $@ 

