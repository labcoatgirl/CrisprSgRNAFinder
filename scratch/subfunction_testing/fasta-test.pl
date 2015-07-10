

$header= ">hg19_dna range=chr3:8790600-8811731 5'pad=0 3'pad=0 strand=+ repeatMasking=none9";

if ($header =~ /.*(chr\w+):(\d+)-(\d+)/g ){print "$1 $2 $3"}