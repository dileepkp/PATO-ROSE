:- module(readable_loop, [readableLoop/1]).

:- use_module(basic).
:- use_module(variable).
:- use_module(loop).
:- use_module(type).

%%	The top level specification
	
readableLoop(LoopURL) :-
	forLoop(LoopURL),
	hasForInit(LoopURL, Init),
	findall(LoopVar,isreadableInit(Init, LoopVar),List),
	% write(LoopURL),nl,nl,
	% write(List),nl,
	hasForTest(LoopURL, Test),
	findall(LoopVar2,isreadableTest(Test, LoopVar2),List2),
	% write(List2),nl,
	\+ subset(List, List2).

%%	isreadableInit

isreadableInit(InitURL, LoopVar) :-
	%% forLoop(LoopURL), hasForInit(LoopURL, InitURL),
	rosereadableInit(InitURL, LoopVar, _).

%% normal style
rosereadableInit(InitURL, LoopVar, LB) :-
	hasChild(InitURL, Temp), 
	hasChild(Temp, AssignOpURL),
	c_is_a(AssignOpURL, 'AssignOp'), !,
	leftOperand(AssignOpURL, VarRefURL),
	variableRef(VarRefURL),
	get_varDecl(VarRefURL, LoopVar).

%%	isreadableTest
	
isreadableTest(TestURL, LoopVar) :- 
	%% forLoop(LoopURL), hasForTest(LoopURL, TestURL),
	rosereadableTest(TestURL, LoopVar).

rosereadableTest(TestURL, LoopVar) :-
	hasChild(TestURL, Temp), 
	hasChild(Temp, VarRefURL),
	variableRef(VarRefURL),
	get_varDecl(VarRefURL, LoopVar).





