%% #!/usr/bin/env swipl
%% swipl --nosignals --quiet run.pl -- is.c.ttl is.c.report.txt

:- set_prolog_flag(verbose, silent).

:- use_module(library(semweb/rdf_db)).
:- use_module(library(semweb/turtle)).

:- ensure_loaded(load).

:- initialization main.

:- rdf_register_prefix(c, 'http://www.semanticweb.org/yzhao30/ontologies/2015/7/c#').

eval :-
	current_prolog_flag(argv, Argv),
	%% [Inputfile, Outputfile|_] = Argv,
	% append(_, [-- | Allargs], Argv),
	append(_, Allargs, Argv),
	[Inputfile, Outputfile|_] = Allargs,
	%% c_onto(Conto),
	%% rdf_load(Conto, [register_namespaces(ture)]),
	%% writeln(Inputfile), 
 open(Inputfile, read, Stream),
 read_string(Stream, "\n", "\r", End1, String1),
 read_string(Stream, "\n", "\r", End2, String2),
 read_string(Stream, "\n", "\r", End3, String3),
 read_string(Stream, "\n", "\r", End4, String4),
 read_string(Stream, "\n", "\r", End5, String5),
 read_string(Stream, "\n", "\r", End6, String6),
 read_string(Stream, "\n", "\r", End7, String7), close(Stream),
 writeln(String1),writeln(String2),writeln(String3),writeln(String4),writeln(String5),writeln(String6),writeln(String7),
  nl,nl,
	analysis(Inputfile, Outputfile).
		
main :- 
	catch(eval, E, (print_message(error, E), fail)), 
	halt(1).
	
main :-
	halt(1).


analysis(InTurtle, OutReport) :-
	% writeln('-----------------------------------'),
	% write('|    '), write(InTurtle), writeln('    |'),
	% writeln('-----------------------------------'),
	atom(OutReport), open(OutReport, write, Out), !,
	%% load knowledge base of input program
	statistics(cputime, T1),
	rdf_load(InTurtle, [format(turtle), register_namespaces(true)]),
	
	% forall( isFunctionDef(FunctionDef), writeln(FunctionDef) ),
	forall( isFunctionDef(FunctionDef), cfg(FunctionDef) ),

	statistics(cputime, T2), % in seconds
	Tdiff is (T2 - T1),
	% writeln('-----------------'),
	% write('#'), writeln(Tdiff), 
	% writeln('-----------------'),
	close(Out).
	% writeln(Tdiff),
	% writeln('------------Done------------------').


fwriteln(Out, Text) :-
	write(Out, Text), nl(Out).


