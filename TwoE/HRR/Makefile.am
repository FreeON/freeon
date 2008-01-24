include ../../make.inc

noinst_LIBRARIES = libHRR.a

libHRR_a_SOURCES = \
  KetHRR1010.f90 \
  KetHRR102.f90 \
  KetHRR103.f90 \
  KetHRR104.f90 \
  KetHRR106.f90 \
  KetHRR1510.f90 \
  KetHRR1515.f90 \
  KetHRR152.f90 \
  KetHRR153.f90 \
  KetHRR154.f90 \
  KetHRR156.f90 \
  KetHRR21.f90 \
  KetHRR22.f90 \
  KetHRR32.f90 \
  KetHRR33.f90 \
  KetHRR42.f90 \
  KetHRR43.f90 \
  KetHRR44.f90 \
  KetHRR62.f90 \
  KetHRR63.f90 \
  KetHRR64.f90 \
  KetHRR66.f90

libHRR_a_DEPENDENCIES = \
  $(foreach i, $(MondoMods_DEPENDENCY_FILES), ../$(i)) \
	../MSC/shellpairstruct.mod \
	../MSC/vscratchb.mod \
	../MSC/getatompairmod.mod

EXTRA_INCLUDES = -I../../MondoMods -I../MSC
