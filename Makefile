# ----- 项目基本信息 -----
NAME 			= WTA_CG
CXX 			= g++
CFLAGS 			= -O2 -pipe -g

# ------ 用于CPLEX的接口
CPLEXLIBDIR   = cplex/lib/

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CLNDIRS   = -I$(CPLEXINCDIR) -L$(CPLEXLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread -ldl -lhighs
CLNFLAGS  = -lcplex -lm -lpthread -ldl -lhighs

CPLEXINCDIR   = cplex/include




# ----- 给出一些内容的路径 -----
CUR_DIR			= $(shell pwd)

# Objecte files
OBJDIR			= OBJ
BINOBJDIR		= ${OBJDIR}/bin
LIBOBJDIR 		= ${OBJDIR}/lib

# Binary file
BINDIR			= bin

# SRC Files and main.cpp
SRCDIR			= src
MAINDIR			= src

# Data Files
DATADIR			= data


# -------- 枚举头文件与OBJ文件 -----
# 头文件
HEADER	 		= 	Solution.h \
					Node.h \
					AlgorithmParameter.h \
					WTA.h \
					SubProblem.h \
					Scene.h \
					Pricing.h \
					print.h \
					Master.h \
					Tree.h \
					BranchAndCut.h \
# OBJ文件
LIBOBJ 			= 	WTA.o \
					SubProblem.o \
					Pricing.o \
					Master.o \
					BranchAndCut.o \
					Scene.o\

# main函数的OBJ文件
MAINOBJ 		= 	mainWTA.o

# ------ 通过.h 与 .o 文件获得【包含路径的】.cpp文件 与 .h文件 .o文件 与 二进制文件
#.cpp 文件
LIBSRC			= $(addprefix $(SRCDIR)/,$(LIBOBJ:.o=.cpp))
MAINSRC     	= $(addprefix $(MAINDIR)/,$(MAINOBJ:.o=.cpp))

# .h文件
LIBSRCHEADER 	= $(addprefix $(SRCDIR)/,$(HEADER))

# .o文件
MAINOBJFILES 	= $(addprefix $(BINOBJDIR)/,$(MAINOBJ))
LIBOBJFILES 	= $(addprefix $(LIBOBJDIR)/,$(LIBOBJ))

# 二进制文件
BINFILE 		= $(BINDIR)/$(NAME)


# 全部内容：生成所需路径、生成库文件与二进制文件
# 我还不了解库的含义，目前尚未加入库文件
all:$(LIBDIR) $(BINDIR) $(OBJDIR) $(LIBOBJDIR) $(BINOBJDIR) $(BINFILE)
	@-$(MAKE) $(LIBDIR) $(BINDIR) $(OBJDIR) $(LIBOBJDIR) $(BINOBJDIR)
	
# 生成二进制文件：依赖于全部的OBJ文件与.h文件
$(BINFILE):$(LIBOBJFILES) $(MAINOBJFILES) $(LIBSRCHEADER)
	g++ $(LIBOBJFILES) $(LPIOBJFILE) $(MAINOBJFILES) -o $@ $(CLNFLAGS)

# 确保所有的文件路径都已经存在
$(LIBDIR):
	@-mkdir -p $(LIBDIR)
$(BINDIR):
	@-mkdir -p $(BINDIR)
$(OBJDIR):
	@-mkdir -p $(OBJDIR)
$(LIBOBJDIR):
	@-mkdir -p $(LIBOBJDIR)
$(BINOBJDIR):
	@-mkdir -p $(BINOBJDIR)

# 编译.cpp文件为.o文件，包括 main 与其他 src
$(LIBOBJDIR)/%.o:$(SRCDIR)/%.cpp 
	$(CXX) -c $(CFLAGS)  $< -o $@ $(CLNFLAGS)
$(MAINOBJFILES):$(MAINSRC) 
	$(CXX) -c $(CFLAGS)  $< -o $@ $(CLNFLAGS)

#clean
.PHONY:clean
clean:
	@-rm -rf $(LIBDIR)
	@-rm -rf $(BINDIR)
#这里必须先删除两个子目录，不然会出错
	@-rm -rf $(LIBOBJDIR)
	@-rm -rf $(BINOBJDIR)
	@-rm -rf $(OBJDIR)
	@-rm  -f ./bin/$(BINFILE)




