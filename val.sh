make -j
valgrind --log-file="./log/valgrind_report.log" --tool=memcheck --leak-check=yes --track-origins=yes ./bin/WTA_CG
