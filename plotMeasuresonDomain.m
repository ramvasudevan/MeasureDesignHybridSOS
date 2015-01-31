function C = plotMeasuresOnDomain(argument,W,x,Rx,color)
%disp('Pcontour');

global slice4DAt

hold on;
gridStep=0.01;

if isempty(slice4DAt)
    slice4DAt=zeros(4,1);
end

if ~exist('color', 'var')
    if strcmp(argument, 'W')
        color = 'r';
    elseif strcmp(argument, 'V')
        color = 'b';
    end
    
end

%disp('p');
linewidth=4;
linespec = color;
% linespec=['-' color];

switch size(x,1)    %Plot Based on the dimention
    case 1
        temp = [];
        if numel(Rx) ~= 2
            disp('ERROR: WRONG NUMBER OF RX');
        end
        x1 = Rx(1):0.01:Rx(2);
        for i=1:numel(x1)
            temp(i) = subs(W, x, x1(i));
        end
        plot(x1, temp, linespec); hold on;
        plot([Rx(1) -0.5 -0.5 0.5 0.5 Rx(2)], [0 0 1 1 0 0], 'r');
        hold off;
    case 2
        
        if isa(W, 'sym')
            problem = W;
        else
            x = sym('x',[2 1]);
            problem=spot2syms(W, x);
            if numel(Rx) ~= 4
                disp('ERROR: WRONG NUMBER OF RX');
            end            
        end
        domain = [-1 1 -1 1] ;
        [y1,y2] = meshgrid(Rx(1,1):gridStep:Rx(1,2),Rx(2,1):gridStep:Rx(2,2));
        
        hold on;
        if strcmp(argument, 'W')
            pgrid = double(subs(problem,{'x1','x2'},{y1,y2}));
            [c,h] = contour(y1,y2,pgrid,[1,1],linespec);
        elseif strcmp(argument, 'V')
            pgrid = double(subs(problem,{'x1','x2','t'},{y1,y2,0}));
            [c,h] = contour(y1,y2,pgrid,[0,0],linespec);
        end
        set(h,'linewidth',2);
        
        xlabel('x1');
        ylabel('x2');
        axis (domain);
        
        y1V = y1(:);
        y2V = y2(:);
        val = pgrid(:);
        idx = find(val >= 1);
        C.samplePts = [y1V(idx), y2V(idx)];
        C.edges = c;
        
        idxNotIn = find(val >= 0.95 & val < 1);
        C.narrowlyNotIn = [y1V(idxNotIn), y2V(idxNotIn)];
    case 4
        w=sdisplay(W);
        x=sdisplay(x);
        if numel(Rx) == 1
            Rx = [Rx Rx];
        end
        domain = [-Rx(1) Rx(1) -Rx(2) Rx(2)] ;
        %disp('2D Problem');
        [y1,y2] = meshgrid(-Rx(1):gridStep:Rx(1),-Rx(2):gridStep:Rx(2));
        problem = w{1};
        problem = strrep(problem,'x(1)','x1');
        problem = strrep(problem,'x(2)','x2');
        problem = strrep(problem,'x(3)','0');
        problem = strrep(problem,'x(4)','0');
        problem =  sym(problem);
        
        hold on;
        C = [];
        if strcmp(argument, 'W')
            pgrid = double(subs(problem,{'x1','x2'},{y1,y2}));
            contour(y1,y2,pgrid,[1,1],linespec);
        elseif strcmp(argument, 'V')
            pgrid = double(subs(problem,{'x1','x2','t'},{y1,y2,0}));
            contour(y1,y2,pgrid,[0,0],linespec);
        end
        %pgrid1 = pgrid >= 1 ;
        %contourf(y1,y2,pgrid1,[1,1],linespec);
        % use contourf?
        C.y1 = y1;
        C.y2 = y2;
        C.pgrid = pgrid;
        xlabel('x1');
        ylabel('x2');
        axis (domain);
        
        %         w=sdisplay(W);
        %         x=sdisplay(x);
        %         if numel(Rx) == 1
        %             Rx = [Rx Rx];
        %         end
        %         domain = [-Rx(1) Rx(1) -Rx(2) Rx(2)] ;
        %         %disp('2D Problem');
        %         [y1,y2,y3,y4] = ndgrid(-Rx(1):gridStep:Rx(1),-Rx(2):gridStep:Rx(2),-Rx(3):gridStep:Rx(3),-Rx(4):gridStep:Rx(4) );
        %         problem = w{1};
        %         problem = strrep(problem,'x(1)','x1');
        %         problem = strrep(problem,'x(2)','x2');
        %         problem = strrep(problem,'x(3)','x3');
        %         problem = strrep(problem,'x(4)','x4');
        %         problem =  sym(problem);
        %
        %         hold on;
        %         C = [];
        %         if strcmp(argument, 'W')
        %             pgrid = double(subs(problem,{'x1','x2','x3','x4'},{y1,y2,y3,y4}));
        %             contour(y1,y2,pgrid,[1,1],linespec);
        %         elseif strcmp(argument, 'V')
        %             pgrid = double(subs(problem,{'x1','x2','x3','x4','t'},{y1,y2,y3,y4,0}));
        %             contour(y1,y2,pgrid,[0,0],linespec);
        %         end
        %         %pgrid1 = pgrid >= 1 ;
        %         %contourf(y1,y2,pgrid1,[1,1],linespec);
        %         % use contourf?
        %         C.y1 = y1;
        %         C.y2 = y2;
        %         C.pgrid = pgrid;
        %         xlabel('x1');
        %         ylabel('x2');
        %         axis (domain);
end

