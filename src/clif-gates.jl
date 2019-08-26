

clifX = [1 0 1; 0 1 0; 0 0 1]

clifY = [1 0 1; 0 1 1; 0 0 1]

clifZ = [1 0 0; 0 1 1; 0 0 1]

clifHpart = [0 1 0; 1 0 0; 0 0 1]
# clifH also needs to update s as s + z*x

clifSpart = [1 0 0; 1 1 0; 0 0 1]
# clifS also needs to update s as s + z*x


#clifCNOT sends
# ZI → ZI
# IZ → ZZ
# XI → XX
# IX → IX
clifCNOTpartZ = [1 1; 0 1]
clifCNOTpartX = [1 0; 1 1]

#clifCZ sends
# XI → XX
# IX → XX
clifCZpartX = [1 1; 1 1]
