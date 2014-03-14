B=.
LFLAGS= -lm
CFLAGS= -O2
PROGS = get_info
INSTALL	:  clean	$(PROGS)
$(PROGS): 
	-$(CC) $(CFLAGS) $(@F).c $(LFLAGS) -o $@
	@echo $(@F) installed in $B
clean::
	rm -f $(PROGS)
