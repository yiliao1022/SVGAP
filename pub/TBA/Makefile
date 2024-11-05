CC = gcc
CFLAGS = -Wall -Wextra -Werror
CFLAGS += -O0
#CFLAGS += -ggdb3
#CFLAGS += -pg

PROGS = all_bz lav2maf pair2tb single_cov2 blastzWrapper maf_project maf2lav maf2fasta maf_order mafFind get_standard_headers maf_checkThread tba roast multic multiz maf_sort get_covered
  
ARCH ?= $(shell arch)
INSTALL = install -c
INSTALLDIR = /depot/apps/$(ARCH)/bin
DATE=$(shell date +"%m%d%y")

all : $(PROGS)

pair2tb : pair2tb.c mz_scores.c mz_scores.h maf.c maf.h util.c util.h seq.c seq.h nib.c nib.h charvec.c charvec.h multi_util.c multi_util.h maftop2tb.h maftop2tb.c
	$(CC) $(CFLAGS) pair2tb.c mz_scores.c multi_util.c maf.c util.c seq.c nib.c charvec.c maftop2tb.c -o pair2tb

single_cov2 : single_cov2.c mz_scores.c maf.c maf.h util.c util.h mz_scores.h maf.h seq.c seq.h nib.c nib.h charvec.c charvec.h multi_util.c multi_util.h
	$(CC) $(CFLAGS) single_cov2.c mz_scores.c multi_util.c maf.c util.c seq.c nib.c charvec.c -o single_cov2

lav2maf : lav2maf.c mz_scores.c maf.c maf.h util.c util.h mz_scores.h multi_util.c multi_util.h seq.c seq.h nib.c nib.h charvec.c charvec.h 
	$(CC) $(CFLAGS) lav2maf.c mz_scores.c maf.c util.c seq.c nib.c charvec.c multi_util.c -o lav2maf

all_bz : all_bz.c util.c util.h
	$(CC) $(CFLAGS) all_bz.c util.c -o all_bz

blastzWrapper : blastzWrapper.c nib.c nib.h charvec.c charvec.h multi_util.c multi_util.h seq.c seq.h util.c util.h maf.h maf.c mz_scores.h mz_scores.c
	$(CC) $(CFLAGS) blastzWrapper.c nib.c charvec.c multi_util.c seq.c util.c maf.c mz_scores.c -o blastzWrapper

maf_project : maf_project.c util.c util.h maf.c maf.h multi_util.c multi_util.h mz_scores.c mz_scores.h seq.c seq.h charvec.c charvec.h nib.c nib.h maf_order.h maf_order.c
	$(CC) $(CFLAGS) maf_project.c util.c maf.c multi_util.c mz_scores.c seq.c charvec.c nib.c maf_order.c -o maf_project

maf2lav : maf2lav.c util.h util.c multi_util.h multi_util.c seq.h seq.c nib.h nib.c charvec.h charvec.c maf.h maf.c mz_scores.h mz_scores.c
	$(CC) $(CFLAGS) maf2lav.c util.c multi_util.c seq.c nib.c charvec.c maf.c mz_scores.c -o maf2lav

maf2fasta : maf2fasta.c util.h util.c maf.h maf.c multi_util.h multi_util.c seq.h seq.c nib.h nib.c charvec.h charvec.c mz_scores.h mz_scores.c
	$(CC) $(CFLAGS) maf2fasta.c util.c maf.c multi_util.c seq.c nib.c charvec.c mz_scores.c -o maf2fasta

maf_order : util.h util.c maf.h maf.c multi_util.h multi_util.c maf_order.h maf_order.c seq.h seq.c nib.h nib.c charvec.h charvec.c mz_scores.h mz_scores.c maf_order_main.c
	$(CC) $(CFLAGS) util.c maf.c multi_util.c maf_order.c seq.c nib.c charvec.c mz_scores.c maf_order_main.c -o maf_order
	
mafFind : util.h util.c maf.h maf.c multi_util.h multi_util.c mz_scores.h mz_scores.c mafFind.c
	$(CC) $(CFLAGS) util.c maf.c multi_util.c mz_scores.c mafFind.c -o mafFind

get_standard_headers : util.h util.c seq.h seq.c nib.h nib.c charvec.h charvec.c get_standard_headers.c
	$(CC) $(CFLAGS) util.c seq.c nib.c charvec.c get_standard_headers.c -o get_standard_headers

maf_checkThread : util.h util.c multi_util.h multi_util.c seq.c seq.h nib.h nib.c charvec.h charvec.c mz_scores.h mz_scores.c maf.h maf.c maf_checkThread.c
	$(CC) $(CFLAGS) util.c multi_util.c seq.c nib.c charvec.c mz_scores.c maf.c maf_checkThread.c -o maf_checkThread

roast : util.h util.c multi_util.h multi_util.c maf.h maf.c nib.h nib.c seq.h seq.c charvec.h charvec.c mz_scores.h mz_scores.c speciesTree.h speciesTree.c auto_mz.c
	$(CC) $(CFLAGS) util.c multi_util.c maf.c nib.c seq.c charvec.c mz_scores.c auto_mz.c speciesTree.c -o roast

tba : util.h util.c align_util.h align_util.c multi_util.h multi_util.c maf.h maf.c mz_scores.h mz_scores.c speciesTree.h speciesTree.c nib.h nib.c seq.h seq.c charvec.h charvec.c mz_preyama.h mz_preyama.c mz_yama.h mz_yama.c tba.c
	$(CC) $(CFLAGS) util.c multi_util.c maf.c nib.c seq.c charvec.c mz_scores.c speciesTree.c mz_preyama.c mz_yama.c tba.c align_util.c -o tba

multic : util.h util.c maf.h maf.c multi_util.h multi_util.c mz_scores.h mz_scores.c mz_preyama.h mz_preyama.c mz_yama.h mz_yama.c multic.c seq.h seq.c nib.h nib.c charvec.h charvec.c align_util.h align_util.c
	$(CC) $(CFLAGS) util.c multi_util.c maf.c mz_scores.c nib.c seq.c charvec.c mz_yama.c mz_preyama.c align_util.c multic.c -o multic

multiz : util.h util.c multi_util.h multi_util.c maf.h maf.c mz_scores.h mz_scores.c speciesTree.h speciesTree.c nib.h nib.c seq.h seq.c charvec.h charvec.c multiz.c mz_preyama.h mz_preyama.c mz_yama.h mz_yama.c align_util.h align_util.c
	$(CC) $(CFLAGS) util.c multi_util.c maf.c mz_scores.c speciesTree.c nib.c seq.c charvec.c mz_preyama.c mz_yama.c multiz.c align_util.c -o multiz

maf_sort : util.h util.c maf.h maf.c multi_util.h multi_util.c maf_sort.h maf_sort.c maf_sort_main.c seq.h seq.c nib.h nib.c charvec.h charvec.c mz_scores.h mz_scores.c
	$(CC) $(CFLAGS) util.c maf.c multi_util.c maf_sort.c maf_sort_main.c seq.c nib.c charvec.c mz_scores.c -o maf_sort

get_covered : util.h util.c multi_util.h multi_util.c maf.h maf.c mz_scores.h mz_scores.c get_covered.c seq.h seq.c nib.h nib.c charvec.h charvec.c
	$(CC) $(CFLAGS) util.c multi_util.c maf.c mz_scores.c get_covered.c seq.c nib.c charvec.c -o get_covered

tarfile :
	@mkdir multiz-tba.$(DATE)
	@cp README* Makefile *.c *.h multiz-tba.$(DATE)
	@tar cfz multiz-tba.$(DATE).tar.gz multiz-tba.$(DATE)
	@rm -rf multiz-tba.$(DATE)
	@echo Built tarfile multiz-tba.$(DATE).tar.gz

install :
	mkdir -p $(INSTALLDIR); set -x; for f in $(PROGS); do $(INSTALL) $$f $(INSTALLDIR) || exit 1; done

clean:
	$(RM) $(PROGS) multiz-tba.$(DATE).tar.gz
