# The buildmaster files.
BUILDMASTER_DIR = buildmaster
BUILDMASTER_FILES = \
  $(BUILDMASTER_DIR)/README \
  $(BUILDMASTER_DIR)/master.cfg \
  $(BUILDMASTER_DIR)/buildbot.tac

# Copy buildmaster files to buildmaster.
buildmaster-local :
	rsync -av $(BUILDMASTER_FILES) /home/buildmaster/

buildmaster :
	@echo "FIXME"

clean :
	rm -f *~
