all:
	cd MakeTables ; $(MAKE)
	cd PONXMults  ; $(MAKE)
	cd PONX       ; $(MAKE)

clean:
	cd MakeTables ; $(MAKE) -i clean
	cd PONXMults  ; $(MAKE) -i clean
	cd PONX       ; $(MAKE) -i clean
	cd Drivers    ; rm -f *.binary

