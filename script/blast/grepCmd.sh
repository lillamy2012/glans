for i in $(grep -v "#" out.out | awk '{  print  $2  }'); do
    grep $i ../data/blast/merged.gtf | awk '{print $1,$12,$4,$5}' | awk 'NR == 1 {scf = $1; min = $3; max = $4;tr=$2;}
    NR > 1 && $3 < min {scf = $1; min = $3}
    NR > 1 && $4 > max {scf = $1; max = $4}
    {
    split(scf,a,"_");
    new=sprintf("%05d", a[2])
    }
    END{print "http://localhost:60151/goto?locus=scaffold"new":"min"-"max }'

done



