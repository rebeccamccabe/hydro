function [] = makeVideo(fname, dt, movieVector)

t=clock();
timestamp = sprintf('_%4d_%02d_%02d_%02d_%02d',t(1),t(2),t(3),t(4),t(5));

myWriter = VideoWriter([fname timestamp],'MPEG-4');
myWriter.FrameRate = 1/dt;
open(myWriter);
warning('off')
writeVideo(myWriter,movieVector);
warning('on')
close(myWriter);

end