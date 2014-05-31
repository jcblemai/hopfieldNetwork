###
set terminal qt size 1000,1000 font "Cantarell,10"

set ylabel "Mean Retrival Error"
set xlabel "Flip Ratio"

set title "Mean Retrival Error as a Funtion of the Flip Ratio\n {/*0.5 N=200, P=5, 800 network realisation per flip ratio}"

# plot '.dat' index 0 using 1:2 w l t "" lw 4, "" using 1:2:3 with yerrorbars t "Standart Error of Mean" lw 2
#replot 'data500000step' index 0 using 1:4 with line t "STDev Pressure 500 000 steps" lw 4, "" using 1:4:5 with yerrorbars t "" lw 2


print "press enter to continue"
pause -1

# set ylabel "Mean Retrival Error"
# set xlabel "Pattern Stored"
# 
# set title "Mean Retrival Error as a Function of the number of Pattern Stored"
# plot 'CapacityEstimation.dat' index 0 using 1:2 w l t "" lw 4, "" using 1:2:3 with yerrorbars t "Standart Error of Mean" lw 2

set title "Transition Time as a Function of Tau, for Different Lambda"
set ylabel "Transition time [step]"
set xlabel "Tau"
plot 'transitionTimeGoodData.dat' index 0 using 2:3 w l t "lambda = 1.3" lw 4
replot 'transitionTimeGoodData.dat' index 1 using 2:3 w l t "lambda = 1.7" lw 4
replot 'transitionTimeGoodData.dat' index 2 using 2:3 w l t "lambda = 2.2" lw 4
replot 'transitionTimeGoodData.dat' index 3 using 2:3 w l t "lambda = 2.5" lw 4

print "press enter to continue"
pause -1

set title "Overlap for all pattern as a function of time"
set ylabel "Overlap [Normalized]"
set xlabel "time"
set xrange [0:100]
set yrange [0:1.3]
plot 'Sequencial.dat' index 0 using :2 w l t "Pattern 0" lw 2
replot 'Sequencial.dat' index 0 using :4 w l t "Pattern 1" lw 2
replot 'Sequencial.dat' index 0 using :6 w l t "Pattern 2" lw 2
replot 'Sequencial.dat' index 0 using :8 w l t "Pattern 3" lw 2
replot 'Sequencial.dat' index 0 using :10 w l t "Pattern 4" lw 2

print "press enter to finish"
pause -1