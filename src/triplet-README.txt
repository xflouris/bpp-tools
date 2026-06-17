
bpp4.8 --simulate MCcoal-3s-msc.ctl
bpp4.8 --simulate MCcoal-3s-inflow.ctl

bpp-tools --triplet --msa mydata-3s-msc.txt --map imap-3s.txt

bpp-tools --triplet --msa mydata-4s.txt --map imap-4s.txt
bpp-tools --debug --triplet --msa mydata-4s.txt --map imap-4s.txt

bpp-tools --debug --dstat "A,B,C,D" --msa mydata-4s.txt --map imap-4s.txt


Some notes

inside triplet.c

   int transform_triplet = 1;

this can be changed to 0 (no transform), 1 (logit), and 2 (logarithm).

inside triplet-models.c, 

   opt_debug = 1;

this flags allow intermediate results to be printed out.
