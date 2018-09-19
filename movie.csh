set i = 0
set sz = `./zeroseq $i 5`
while( -e $sz.ps )
	echo $sz.ps $sz.jpg
	ps2pdf -dEPSCrop $sz.ps
	convert $sz.ps $sz.jpg

	set i = `echo $i 1 + p | dc`
	set sz = `./zeroseq $i 5`
end
