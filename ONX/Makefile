all:
#	cd MakeTables ; $(MAKE)
#	cd PONXMults  ; $(MAKE)
#	cd PONX       ; $(MAKE)
	cd NONX       ; $(MAKE)
#
release:clean
	$(MAKE) -C NONX release
	rm -rf MakeTables PONXMults PONX GMEONX GSONX MEONX \
               Drivers GDrivers VRR VRRm ONX SONX MASONX
#
clean:
#	cd MakeTables ; $(MAKE) -i clean
#	cd PONXMults  ; $(MAKE) -i clean
#	cd PONX       ; $(MAKE) -i clean
	cd NONX       ; $(MAKE) -i clean
#	cd Drivers    ; rm -f *.binary
