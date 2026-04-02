#Binary data analysis
for  ((i=1;i<=1;i++))
    	do
		R CMD BATCH --no-save gweis_p.R
		R CMD BATCH --no-save gweis_r.R
		R CMD BATCH --no-save gcim_gweis_p.R
		R CMD BATCH --no-save gcim_gweis_r.R
	done
