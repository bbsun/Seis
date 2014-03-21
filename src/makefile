include OBJX
all:
	@set -e; \
	./mk
	rm makefile_3
.PHONY:clean
clean:
	find . -name "*.o" |xargs rm -f
	find . -name "*.d" |xargs rm -f
	find . -name "*~"  |xargs rm -f
	find . -name "\#*" |xargs rm -f
	find . -name "*\#" |xargs rm -f
	rm run -f	
