include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/RemoveAll
#
inpts:
	./mkinpts
#
run:
	mondo c20.inp  c20.out  c20.log
	mondo c30.inp  c30.out  c30.log    
	mondo c50.inp  c50.out  c50.log    
	mondo c60.inp  c60.out  c60.log    
	mondo c80.inp  c80.out  c80.log
	mondo c180.inp c180.out c180.log    
	mondo c240.inp c240.out c240.log     
#	mondo c320.inp c320.out c320.log   
#	mondo c540.inp c540.out c540.log    
#
lsf:
#
pbs:
	sub_pbs c20
	sub_pbs c30
	sub_pbs c50
	sub_pbs c60
	sub_pbs c80
	sub_pbs c180
	sub_pbs c240
#	sub_pbs c320
#	sub_pbs c540
#
purge:	clean
	rm -f *.inp
	rm -f *.pbs
#
clean:
	rm -f $(REMOVEALL) 
	rm -f *.out *.log *.o*
	rm -f *.inf.ascii
	rm -rf  core core.* core*
	rm -f PI*
	rm -f *.usertime.* *.pcsamp.* *.*_Info64 *_Info64 *.aggregate.* 
	ln -s /dev/null core
