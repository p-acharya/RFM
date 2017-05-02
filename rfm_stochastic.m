% This function performs a stochastic simulation of a TASEP like model
% with parallel update mode. Results are then avaraged to estimate RFM results.
%
%  Usage: [ occupancies, delays ]  = rfm_stochastic( lambda, time_step, sim_time );
%
%  Where:  lambda - an n+1 vector of lambda_i (n is the number of sites).
%          time_step - simulation tick. 
%          sim_time - simulation duration (total of sim_time/time_step ticks).
%
%          occupancies - average occupancies per site.
%          delays - a vector of particle's delays. Vector length is the total
%                   number of particles that traverse the chain of n sites.
%
%
% Yoram Zarai, 10/27/14

% -----------------------------------------------------------------------------------------------------

function [ occupancies, delays ]  = rfm_stochastic( lambda, time_step, sim_time );
  
n = length( lambda ) - 1; % number of nodes
tlambda = [ lambda, 1 ]; % last entry value is not actually used (this is done to 
                         % save an if statement in the loop)

% 'o' - occupancy, 'id' - ID, 'et' - entry time
% node_state.o(2:end-1) contains the nodes' state variables (x_1 to x_n)
% node_state.o( 1 ) always contains 1, node_state.o( end ) always contains 0
% we refer to node_state.x(1) as node 0 and node_state.x(end) as sink node. Both are
% not part of the RFM/TASEP n nodes chain.
node_state = struct( 'o', [1,zeros(1,n),0], 'id', [1 zeros(1,n+1) ], 'et', zeros(1,n+2) );

occu = zeros( 1, length( node_state.o ) );
time_to_hop = [0 inf*ones(1,n+1)]; % so node 0 will hop at first tick of simulation (so that
                                 % node 1 will contain a particle after the first tick)

id = 2; % next particle ID to be insert into node 0 (first ID is at node 0)
dlys = zeros( 1, round( sim_time/time_step ) );  % measuring delay per particle ID
last_id = 0;

% particle flow simulation
% -------------------------
fprintf( 1, '%s: Counting time....', mfilename ); tic;
%for time = 0 : time_step : sim_time
time = 0;
while( time < sim_time )
  for node = n+1 : -1 : 1  % start from node n up to node 0
    if( time == time_to_hop( node ) )
      
      if( ( node_state.o(node) == 1 ) && (node_state.o(node+1) == 0 ) )
        % particle moves from node to node+1
        node_state.o(node+1) = 1; node_state.id(node+1)=node_state.id(node); node_state.et(node+1)=node_state.et(node); 
        node_state.o(node) = 0; node_state.id(node)=0; node_state.et(node)=0; % last two just for readability
        
        % set time to hop of node that acquired a particle based on the corresponding lambda.
        % node(i) hops with rate tlambda(i).
        time_to_hop( node + 1 ) = max( ( time+time_step ), ( time + exprnd( 1 / tlambda( node + 1 ) ) ) );
        time_to_hop( node + 1 ) = round( time_to_hop( node + 1 ) / time_step ) * time_step; 
      else        
        % particle did not move from node, setting new time to hop
        time_to_hop( node ) = max( ( time+time_step ), ( time + exprnd( 1 / tlambda( node ) ) ) );
        time_to_hop( node ) = round( time_to_hop( node ) / time_step ) * time_step; 
      end;
      
      % node 0 should always have the next particle
      if( node_state.o(1) == 0 ) 
        time_to_hop( 1 ) = max( ( time+time_step ), ( time + exprnd( 1 / tlambda( 1 ) ) ) );
        time_to_hop( 1 ) = round( time_to_hop( 1 ) / time_step ) * time_step; 

        node_state.o( 1 ) = 1; node_state.id(1) = id; id = id+1;
      end;
      
      % node sink should always be empty. If it has a particle, log its statistics and empty it
      if( node_state.o( end ) == 1 )
        node_state.o( end ) = 0;
        last_id = node_state.id( end ); 
        dlys( last_id  ) = time - node_state.et( end );
      end;

      node_state.et(1) = time_to_hop(1); % we don't count for node 0 delay 
    end; % time == time_to_hop
  end; % node loop
  
  occu = occu + node_state.o;
  
  % jump to next time that contains a hop (probably useful for small lambda)
  %next_t = min( sim_time, min( time_to_hop( time_to_hop > time ) ) ) - time_step; % the loop adds time_step
  %occu = occu + ((next_t - time )/time_step) * node_state.o;
  %time = round( ( next_t + time_step ) / time_step ) * time_step;
  
  time = round( ( time + time_step ) / time_step ) * time_step;
end;
fprintf( 1, 'Done. Run time = %f sec\n', toc );

% averaging 
occupancies = occu( 2 : end - 1 ) / ( sim_time / time_step );
delays = dlys( 1 : last_id );

