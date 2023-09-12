function a = angle(x,y)
len = length(x(:));
for i = 1:len
          if atan2(y(i),x(i))>=0
                a(i) = (180/pi) * atan2(y(i),x(i));
          else
                 a(i) = (180/pi) * (atan2(y(i),x(i))+2*pi);
           end
end
end