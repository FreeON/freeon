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
	rm -rf PONXMults PONX GMEONX GSONX MEONX \
               Drivers GDrivers VRRm ONX SONX MASONX
#
purge:	clean
	make -C NONX purge
	rm -rf Drivers/*.binary
#
clean:
	cd VRR        ; $(MAKE) -i clean
	cd NONX       ; $(MAKE) -i clean
#	cd MakeTables ; $(MAKE) -i clean
#	cd Drivers    ; rm -f *.binary
