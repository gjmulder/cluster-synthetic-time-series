#!/bin/sh

tail -n 100000 -F nohup.out | mawk -W interactive '
/Number of combination targets/ {
	combs=$(NF-1);
}

/^GA/ {
	printf(".");
}

/^.1. .targets/ {
	count +=1;
	str = "date -Isec";
	str | getline date;
	close(str);
	$1 = "";
	printf("\n%s %40s, combinations: %3d/%3d, progress:%3d%%\n", date, $0, count, combs, 100*count/combs);
}
'
