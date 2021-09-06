function [] = makeVideo(fname, dt, movieVector)

myWriter = VideoWriter(fname,'MPEG-4');
myWriter.FrameRate = 1/dt;
open(myWriter);
writeVideo(myWriter,movieVector);
close(myWriter);

end