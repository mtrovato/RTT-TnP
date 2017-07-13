#!/bin/bash

echo "Building CPlot..."
g++ -O3 -g -c -Wall -fPIC -o CPlot.o CPlot.cc  `root-config --cflags` `root-config --libs` -I $ROOFITSYS/include -L $ROOFITSYS/lib -l RooFit -l RooFitCore
g++ -O3 -shared -Wl,-soname,libCPlot.so -o libCPlot.so CPlot.o

echo "Building KStyle..."
g++ -O3 -g -c -Wall -fPIC -o KStyle.o KStyle.cc `root-config --cflags` `root-config --libs`
g++ -O3 -shared -Wl,-soname,libKStyle.so -o libKStyle.so KStyle.o

# echo "Building RooCMSShape..."
# g++ -O3 -g -c -Wall -fPIC -o RooCMSShape.o RooCMSShape.cc `root-config --cflags` `root-config --libs` -I $ROOFITSYS/include -L $ROOFITSYS/lib -l RooFit -l RooFitCore
# g++ -O3 -shared -Wl,-soname,libRooCMSShape.so -o libRooCMSShape.so RooCMSShape.o

# echo "Building RooVoigtianShape..."
# g++ -O3 -g -c -Wall -fPIC -o RooVoigtianShape.o RooVoigtianShape.cc `root-config --cflags` `root-config --libs` -I $ROOFITSYS/include -L $ROOFITSYS/lib -l RooFit -l RooFitCore
# g++ -O3 -shared -Wl,-soname,libRooVoigtianShape.so -o libRooVoigtianShape.so RooVoigtianShape.o

echo "Building CEffUser1D..."
g++ -O3 -g -c -Wall -fPIC -o CEffUser1D.o CEffUser1D.cc `root-config --cflags` `root-config --libs`
g++ -O3 -shared -Wl,-soname,libCEffUser1D.so -o libCEffUser1D.so CEffUser1D.o

# echo "Building CEffUser2D..."
# g++ -O3 -g -c -Wall -fPIC -o CEffUser2D.o CEffUser2D.cc `root-config --cflags` `root-config --libs`
# g++ -O3 -shared -Wl,-soname,libCEffUser2D.so -o libCEffUser2D.so CEffUser2D.o

echo "Building EffRTTFitter..."
g++ -O3 -g -c -Wall -fPIC -o EffRTTFitter.o EffRTTFitter.cc `root-config --cflags` `root-config --libs` -l Minuit -I $ROOFITSYS/include -L $ROOFITSYS/lib -l RooFit -l RooFitCore -L. -l CPlot -l KStyle -l CEffUser1D
#-l RooCMSShape -l RooVoigtianShape -l CEffUser1D -l CEffUser2D
g++ -O3 -shared -Wl,-soname,libEffRTTFitter.so -o libEffRTTFitter.so EffRTTFitter.o

echo "Building effRTTFit..."
g++ -O3 -g -o effRTTFit effRTTFit.cc `root-config --cflags` `root-config --libs` -l Minuit -I $ROOFITSYS/include -L $ROOFITSYS/lib -l RooFit -l RooFitCore -L. -l EffRTTFitter -l CPlot -l KStyle -l CEffUser1D
#-l RooCMSShape -l RooVoigtianShape -l CEffUser1D -l CEffUser2D
