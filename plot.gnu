###
set terminal qt size 1000,1000 font "Cantarell,10"

set ylabel "Mean Retrival Error"
set xlabel "Flip Ratio"

set title "Mean Retrival Error as a Funtion of the Flip Ratio\n {/*0.5 N=200, P=5, 800 network realisation per flip ratio}"

# plot '.dat' index 0 using 1:2 w l t "" lw 4, "" using 1:2:3 with yerrorbars t "Standart Error of Mean" lw 2
#replot 'data500000step' index 0 using 1:4 with line t "STDev Pressure 500 000 steps" lw 4, "" using 1:4:5 with yerrorbars t "" lw 2



print "press enter to continue"
pause -1

set ylabel "Mean Retrival Error"
set xlabel "Pattern Stored"

set title "Mean Retrival Error as a Function of the number of Pattern Stored"
plot 'CapacityEstimation.dat' index 0 using 1:2 w l t "" lw 4, "" using 1:2:3 with yerrorbars t "Standart Error of Mean" lw 2


print "press enter to finish"
pause -1