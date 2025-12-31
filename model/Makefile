#Makefile for building and running the jjm model on a Linux or Mac OSX
EXEC1=jjm
EXEC2=jjms
# TPL=$(EXEC1).tpl

all: $(EXEC1) $(EXEC2)

$(EXEC1): $(EXEC1).tpl
	admb -f $(EXEC1)

$(EXEC2): $(EXEC2).tpl
	admb -f $(EXEC2)

clean:
	rm -f $(EXEC1) $(EXEC2) *.eva *.log *.rpt *.htp *.cor *.par *.r* *.b* *.p* *.obj
	rm -f *.rep *.bar *.psv *.std $(EXEC2).cpp $(EXEC1).cpp admodel.* variance *.o *.ecm *.mc2 *.mcm
