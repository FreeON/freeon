include $(srcdir)/../validate.inc

bin_SCRIPTS = \
  MDO2xyz.py \
  get_energies.py

VALIDATE_INPUTS = \
  MD_C2F4_Symplectic.inp \
  MD_C2F4_TRBO_ACT_2.inp \
  MD_C2F4_TRBO_ACT.inp \
  MD_C2F4_TRBO.inp \
  MD-F2-Symplectic.inp \
  MD-F2-TRBO.inp \
  MD-H2-BO.inp \
  MD-H2.inp \
  MD-H2O.inp \
  MD-H2O-TRBO-BasisSwitch.inp \
  MD-H2O-TRBO-Damp_dt11.inp \
  MD-H2O-TRBO-Damp_dt7.inp \
  MD-H2O-TRBO-Damp_dt9.inp \
  MD-H2O-TRBO.inp \
  MD-H2O-TRBO-no_damping.inp \
  MD-H2-Symplectic.inp \
  MD_H2_TRBO_2.inp \
  MD_H2_TRBO_3.inp \
  MD_H2_TRBO_ACT_2.inp \
  MD_H2_TRBO_ACT.inp \
  MD_H2_TRBO_Damp_dt7.inp \
  MD-H2-TRBO.inp \
  MD_Water10_T050_1SCF.inp

VALIDATE_TESTS = $(foreach i, $(VALIDATE_INPUTS), $(basename $i).test)
VALIDATE_REFERENCES = $(foreach i, $(VALIDATE_INPUTS), $(basename $i).reference)

EXTRA_DIST = \
  $(bin_SCRIPTS) \
	$(VALIDATE_INPUTS) \
  $(VALIDATE_REFERENCES) \
  $(VALIDATE_TESTS)
