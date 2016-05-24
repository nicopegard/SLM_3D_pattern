x = 1:10;
f = 'A.gif';
for n = 1:0.5:5
y = x.^n;
plot(x,y)
drawnow
frame = getframe(1);
im = frame2im(frame);
[A,map] = rgb2ind(im,256); 
	if n == 1;
		imwrite(A,map,f,'LoopCount',Inf,'DelayTime',1);
	else
		imwrite(A,map,f,'WriteMode','append','DelayTime',1);
	end
end