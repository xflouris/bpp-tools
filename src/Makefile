# Copyright (C) 2021-2023 Tomas Flouri
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
# Department of Genetics, Evolution and Environment,
# University College London, Gower Street, London WC1E 6BT, England

# Profiling
#PROFILING=-g -pg
PROFILING=-g

# Compiler warnings
WARN=-Wall -Wsign-compare -Wextra

AVXDEF=-DHAVE_AVX
AVX2DEF=-DHAVE_AVX2

ifndef CC
CC = gcc
endif
CFLAGS = -D_GNU_SOURCE -g -O3 -msse3 $(AVXDEF) $(AVX2DEF) $(WARN)
LINKFLAGS=$(PROFILING)
LIBS=-lm

PROG=bpp-tools

all: $(PROG)

OBJS=bpp-tools.o util.o arch.o hash.o phylip.o maps.o msa.o dstat.o hardware.o \
     explode.o list.o extract.o remove.o

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $+ $(LIBS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~ $(OBJS) gmon.out $(PROG)
