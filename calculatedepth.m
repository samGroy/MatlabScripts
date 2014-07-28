for i=length(surf)
    if surf(i,1:2)==faults(i,1:2);
        faults=surf(i,3))-faults(i,3);
    i+1;
    end
end
