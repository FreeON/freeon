all:
#	cd MakeTables ; $(MAKE)
#	cd PONXMults  ; $(MAKE)
#	cd PONX       ; $(MAKE)
#	cd GDrivers   ; $(MAKE)
	cd VRR        ; $(MAKE) 
	cd NONX       ; $(MAKE)
#	cd MakeTables ; $(MAKE)
#
release:clean
	$(MAKE) -C NONX release
	$(MAKE) -C VRR  release
	rm -rf PONXMults PONX GMEONX GSONX MEONX \
               GDrivers VRRm ONX SONX MASONX
#
purge:	clean
#vw	make -C NONX purge
	cd VRR        ; $(MAKE) -i clean #vw
	cd NONX       ; $(MAKE) -i clean #vw
	cd MakeTables ; $(MAKE) -i clean #vw
	cd Drivers    ; rm -f *.binary   #vw
#vw	rm -rf Drivers/*.binary
#
clean:
#vw	cd VRR        ; $(MAKE) -i clean
	cd NONX       ; $(MAKE) -i clean
#	cd MakeTables ; $(MAKE) -i clean
#	cd Drivers    ; rm -f *.binary
