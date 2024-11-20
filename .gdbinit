set args -i foo.clq
define reload
	make
	run
end
b main
b mc.hpp:graph
b read_chunk
