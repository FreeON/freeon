include $(MONDO_HOME)/Includes/CompilerSpecifics

inpts:	
	./mkinpts

lsf:
	./submit_lsf_waters

pbs:
	./submit_pbs_waters

clean:
	rm -f $(REMOVEALL90) 
	rm -f *.out *.inp *.lsf *.pbs  *.log *.o*
	rm -f *.inf.ascii
	rm -rf  core core.* core*
	rm -f PI*
	rm -f *.usertime.* *.pcsamp.* *.*_Info64 *_Info64 *.aggregate.* 