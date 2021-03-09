function [walk] = find_walk(edge,strike_threshold,G3,XY3)
  
   
   chain = sort(edge(1,1:2),2);
   node_s = edge(1,1);  
   a = find_possible_edges3(G3,node_s,chain,XY3,strike_threshold);  
   node_t = edge(1,2);
   chain2 = fliplr(chain);
   b = find_possible_edges3(G3,node_t,chain2,XY3,strike_threshold); 
   c = [a fliplr(b)];
   [~,u,~]=unique(c, 'first') ;
   
   walk = c(sort(u));
end

