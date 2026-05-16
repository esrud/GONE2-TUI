MAXLOCI = 2000000
MAXIND  = 2000

CC      = g++
CLUSTERC = /DATA/APPS/gcc/7.2.0/bin/g++

COMMON_FLAGS = -Wall -fopenmp
CFLAGS       = -O3 -DMAXLOCI=$(MAXLOCI) -DMAXIND=$(MAXIND)
FASTFLAGS    = -Ofast

MAC_PATH_OPENMP = /opt/homebrew/Cellar/libomp/19.1.5
MAC_FLAGS = -Wall -Xpreprocessor -fopenmp -lomp \
            -I$(MAC_PATH_OPENMP)/include -L$(MAC_PATH_OPENMP)/lib -std=c++11

OFNAME   = gone2
TUI_OFNAME = gone2-tui
NCURSES_OFNAME = gone-ncurses
TEST_PROGRESS = progress_tui_test
SOURCES  = gone2.cpp lib/*.cpp
NCURSES_LIBS = -lncursesw -lpthread

all: gone

info:
	@echo
	@echo "Compiling GONE2 with MAXLOCI=$(MAXLOCI) and MAXIND=$(MAXIND)"
	@echo
	@echo "Override either limit on the command line:"
	@echo "    make MAXLOCI=10000000 MAXIND=5000 gone"
	@echo
	@echo "Available targets: gone, tui, gone-ncurses, fast, static, arch, debug, cluster, macos, test-progress"
	@echo

gone: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -o $(OFNAME) $(SOURCES)

tui: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -DGONE_TUI -o $(TUI_OFNAME) $(SOURCES)

# Default ncurses build: combo runtime. Reads the input and computes
# d² exactly once, then runs the Ne-estimation GA three times — once
# per algorithm (trunc05_kick, L2, L1+kicks) — switching the
# smoothness penalty and kick schedule at runtime between passes.
# The TUI status line shows the current pass label; once all three
# finish the binary keeps whichever pass produced the lowest
# bin-weighted d² residual, pushes its Ne curve back to the chart,
# and writes _GONE2_Ne / _GONE2_d2 from that winning pass.
#
# Mutation rates use the original empirical defaults
# (frecMut=0.3, efectoMut=0.3). The 3×3 sweep on data/simus (May
# 2026) found these near-optimal at full strength; a static
# frecMut=0.5 and a Rechenberg 1/5 self-adapting efectoMut both
# improved fast-preset RMSE by 5–7 % but the gap closed at full
# strength, so neither replaced the defaults.
# (GA_FREC_MUT / GA_EFECTO_MUT / GA_FREC_MUT_LATERAL remain
# overridable from the command line for future sweeps.)
gone-ncurses: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -DGONE_NCURSES_TUI -DGONE_COMBO \
	    -o $(NCURSES_OFNAME) $(SOURCES) $(NCURSES_LIBS)

fast: info
	$(CC) $(COMMON_FLAGS) $(FASTFLAGS) -o $(OFNAME) $(SOURCES)

static: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -static -o $(OFNAME) $(SOURCES)

arch: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -mavx2 -march=native -o $(OFNAME) $(SOURCES)

debug: info
	$(CC) $(COMMON_FLAGS) $(CFLAGS) -g -o $(OFNAME) $(SOURCES)

cluster: info
	$(CLUSTERC) $(COMMON_FLAGS) $(CFLAGS) -static -o $(OFNAME) $(SOURCES)

macos: info
	$(CC) $(MAC_FLAGS) $(CFLAGS) -o $(OFNAME) $(SOURCES)

test-progress:
	$(CC) $(COMMON_FLAGS) -O0 -g -DGONE_TUI -I. -o $(TEST_PROGRESS) tests/progress_tui_test.cpp lib/progress.cpp lib/tui.cpp
	./$(TEST_PROGRESS)

clean:
	rm -f $(OFNAME) $(TUI_OFNAME) $(NCURSES_OFNAME) $(TEST_PROGRESS) gone2_perf

.PHONY: all info gone tui gone-ncurses fast static arch debug cluster macos test-progress clean
