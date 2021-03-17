$(OBJDIR)/ana.o    :$(MOR_DIR)/code/ana.f    MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/aux.o    :$(MOR_DIR)/code/aux.f    MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/conv.o   :$(MOR_DIR)/code/conv.f   MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/const.o  :$(MOR_DIR)/code/const.f  MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/dump.o   :$(MOR_DIR)/code/dump.f   MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/ei.o     :$(MOR_DIR)/code/ei.f     MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/filter.o :$(MOR_DIR)/code/filter.f MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/lapack.o :$(MOR_DIR)/code/lapack.f MOR LMOR;  $(FC) -c $(FL0) $< -o $@
$(OBJDIR)/legacy.o :$(MOR_DIR)/code/legacy.f MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/pod.o    :$(MOR_DIR)/code/pod.f    MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/qoi.o    :$(MOR_DIR)/code/qoi.f    MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/read.o   :$(MOR_DIR)/code/read.f   MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/riesz.o  :$(MOR_DIR)/code/riesz.f  MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/rk.o     :$(MOR_DIR)/code/rk.f     MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/rom.o    :$(MOR_DIR)/code/rom.f    MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/time.o   :$(MOR_DIR)/code/time.f   MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/tensor.o :$(MOR_DIR)/code/tensor.f MOR LMOR;  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/mpar.o   :$(MOR_DIR)/code/mpar.f   MORDICT MOR LMOR;  $(FC) -c $(FL2) $< -o $@

$(OBJDIR)/unit.o   :$(MOR_DIR)/tests/unit.f  MOR LMOR;  $(FC) -c $(FL2) $< -o $@
