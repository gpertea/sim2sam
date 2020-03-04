# polyester2SAM
conversion of raw polyester FASTA output to SAM in genomic coordinates

#sim2sam Help Page

```
-g -i -o -s -t [-r ]
Arguments:
	g/--gff	annotation from which reads are simulated
	i/--index	Fasta index of the genome from which reads where simulated
	o/--output	name for the output SAM alignment
	r/--rsemi	Base name of the RSEM-prepared reference. Required if the RSEM mode is enabled
	s/--single	single-end reads
	t/--type	Type of the simulator used. Options are [rsem,polyester]
```
