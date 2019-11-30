:- module(cfg, [
	cfg/1,
	cfg/3,
	getEnclosingLoop/2,
	getEnclosingLoopOrSwitch/2,
	getEnclosingFunction/2
]).

:- use_module(basic).
:- use_module(statement).
:- use_module(function).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
%
%	per statement cfg
%	every 
%	+Stmt to statement to analyze
%	-Entry return single entry statement
%	-Exit return single exit (create one if needed)
%
%   @FIXME Goto, Switch
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% cfg for function definition
%	see document.txt for function structure in ROSE   
% cfg(Function) :-
%	functionBodyOf(Function, FunctionBody),functionParameterList(Function, ParameterListURL), !,
% (atom(ParameterListURL) -> parameterProcess(Function)),
 % cfg(FunctionBody, _, _).

cfg(Function) :-
	functionBodyOf(Function, FunctionBody), !,
  parameterProcess(Function,Entry,Exit),
cfg(FunctionBody, _, _, Entry,Exit,Function). 
% ( Entry == none -> cfg(FunctionBody, _, _);
%  cfg(FunctionBody, _, _)
% ).

 parameterProcess(Block, Entry, Exit) :-
	 findall(ChildStmt, hasParameter(Block, _ , ChildStmt), ChildStmts), 
	( % writeln(ChildStmts),
		[] = ChildStmts -> Entry = none, Exit = Block; % empty body
		[First|Rest] = ChildStmts, % otherwise
		cfg(First, FirstEntry, FirstExit),
		% writeln(First),writeln(FirstEntry),
		Entry = FirstEntry, 
		cfg_sq(Rest, FirstExit, Exit)
		% var(Exit) -> Exit = Entry
	).

%%% variable declaration
cfg(VarDecl, Entry, Exit) :-
	(isVariableDecl(VarDecl);isVariable(VarDecl)), !,
	Entry = VarDecl,
	Exit = VarDecl.

%%%  The basic statement (Expression statement)
cfg(Stmt, Entry, Exit) :-
	isBasicStatement(Stmt), !,
	Entry = Stmt,
	Exit = Stmt.

%%% The if-then-else statement
%	create a exit(stmt_id) exit node   
cfg(Stmt, Entry, Exit) :-
	isIfElseStatement(Stmt), !,
	c_has(hasCondition, Stmt, ConditionStmt),
	Entry = ConditionStmt,
	% Entry = Stmt, % dileep
	% Exit = exit(Stmt), % create single blank exit node  dileep
	atomic_list_concat([Stmt, '_exit'], Exit),
	% Exit = Stmt, % dileep
	(
		c_has(hasTrueBody, Stmt, TrueBody) ->
		cfg(TrueBody, TrueEntry, TrueExit),
% write('true '), writeln(TrueEntry),
		setTrueNext(Entry, TrueEntry),
		setNext(TrueExit, Exit) ;
		setTrueNext(Entry, Exit) % if true body is empty
	),
	(
		c_has(hasFalseBody, Stmt, FalseBody) ->
		cfg(FalseBody, FalseEntry, FalseExit),
		setFalseNext(Entry, FalseEntry),
% write('false '), writeln(FalseEntry),
		 setNext(FalseExit, Exit);
		 setFalseNext(Entry, Exit)

	).

%%% The Ternary operator
%	@tbd   

%%% The ForStatement
%	@tbd what if ForInit, ForIncr empty?
% cfg(ForStatement, Entry, Exit) :-
%	isForStatement(ForStatement), !,
%	c_has(hasForInit, ForStatement, Entry), % Entry
%	c_has(hasForTest, ForStatement, ConditionStmt),
%	c_has(hasForIncr, ForStatement, ForIncr),
%	c_has(hasBody, ForStatement, ForBody),
%	Exit = exit(ForStatement), % Exit, blank node 
%	atomic_list_concat([ForStatement, '_exit'], Exit),
%	cfg(ForBody, BodyEntry, BodyExit),
%	setNext(Entry, ConditionStmt),
%	setTrueNext(ConditionStmt, BodyEntry),
%	setFalseNext(ConditionStmt, Exit),
%	setNext(BodyExit, ForIncr),
%	setNext(ForIncr, ConditionStmt).

cfg(ForStatement, Entry, Exit) :-
	isForStatement(ForStatement), !,
	c_has(hasForInit, ForStatement, Entry), % Entry
	c_has(hasForTest, ForStatement, ConditionStmt),
	c_has(hasForIncr, ForStatement, ForIncr),
	c_has(hasBody, ForStatement, ForBody),
%	Exit = exit(ForStatement), % Exit, blank node 
	atomic_list_concat([ForStatement, '_exit'], Exit),
	cfg(ForBody, BodyEntry, BodyExit),
	setNext(Entry, ConditionStmt),
	setTrueNext(ConditionStmt, BodyEntry),
	setFalseNext(ConditionStmt, Exit),
	setNext(BodyExit, ForIncr),
	setNext(ForIncr, ConditionStmt).

%%% The while statement
%	@comment in Rose, while(...); hasBody ExpressionStatement:NullExpr
cfg(Stmt, Entry, Exit) :-
	isWhileStatement(Stmt), !,
	c_has(hasCondition, Stmt, ConditionStmt),
	c_has(hasBody, Stmt, BodyStmt),
	Entry = ConditionStmt, % continue Stmt may point to 
	% Exit = exit(Stmt),% in cfg(BreakStmt), break use the same name to refer to this exit
	atomic_list_concat([Stmt, '_exit'], Exit),
	cfg(BodyStmt, BodyEntry, BodyExit),
	setTrueNext(ConditionStmt, BodyEntry),
	setFalseNext(ConditionStmt, Exit),
	setNext(BodyExit, ConditionStmt),
	setLastStatement(BodyExit,Stmt).

% cfg(Stmt, Entry, Exit) :-
%	isWhileStatement(Stmt), !,
%	c_has(hasCondition, Stmt, ConditionStmt),
%	c_has(hasBody, Stmt, BodyStmt),
%	Entry = Stmt, % continue Stmt may point to 
%	% Exit = exit(Stmt),% in cfg(BreakStmt), break use the same name to refer to this exit
%	atomic_list_concat([Stmt, '_exit'], Exit),
%	cfg(BodyStmt, BodyEntry, BodyExit),
%	setTrueNext(Stmt, BodyEntry),
%	setFalseNext(Stmt, Exit),
%	setNext(BodyExit, Stmt).

%%% The DoWhile ...
cfg(DoWhileStmt, Entry, Exit) :-
	isDoWhileStatement(DoWhileStmt), !,
	c_has(hasBody, DoWhileStmt, BodyStmt),
	c_has(hasCondition, DoWhileStmt, ConditionStmt),
	cfg(BodyStmt, BodyEntry, BodyExit),
	Entry = BodyEntry,
	% Exit = exit(DoWhileStmt),
	atomic_list_concat([DoWhileStmt, '_exit'], Exit),
	setNext(BodyExit, ConditionStmt),
	setTrueNext(ConditionStmt, BodyEntry),
	setFalseNext(ConditionStmt, Exit).

%%% The Switch case statement 
%	@tbd the ROSE gets bug here! cause cycle
cfg(SwitchStmt, Entry, Exit) :-
	isSwitchStatement(SwitchStmt), !,
	Entry = SwitchStmt,
	Exit = SwitchStmt,
	writeln('# SwitchStmt not supported yet').

%%% The Break 
cfg(BreakStmt, Entry, Exit) :-
	isBreakStatement(BreakStmt), !,
	Entry = BreakStmt,
	getEnclosingLoopOrSwitch(BreakStmt, Parent),
	% ParentExit = exit(Parent), % exit(Parent) is a blank exit node, consistent with parent
	atomic_list_concat([Parent, '_exit'], ParentExit),
	setNext(BreakStmt, ParentExit),
	Exit = none. % mark the divergent, will be checked in cfg(Block, ...)
	
%%% The continue 
cfg(ContinueStatement, Entry, Exit) :-
	isContinueStatement(ContinueStatement), !,
	Entry = ContinueStatement,
	getEnclosingLoop(ContinueStatement, Loop),
	(c_has(hasCondition, Loop, NextIterStmt); % while 
		c_has(hasForIncr, Loop, NextIterStmt)), % for loop 
	setNext(ContinueStatement, NextIterStmt),
	Exit = none. % mark the divergent

%%% Goto and labels
%	!not functionaly properly   
cfg(GotoStatement, Entry, Exit) :-
	isGotoStatement(GotoStatement), !,
	Entry = GotoStatement, % !
	Exit = GotoStatement, % !
	writeln('# GotoStatement not supported yet').

%%% Return
%	@tbd should jump to the exit of function 
cfg(ReturnStmt, Entry, Exit) :-
	isReturnStatement(ReturnStmt), !,
	Entry = ReturnStmt,
	getEnclosingFunction(ReturnStmt, FunctionDef),
	% Exit = exit(FunctionDef). dileep
	atomic_list_concat([FunctionDef,'_exit'], Exit).
	  	      
%%% Block contains a sequence of statements 
%	@comment check ChildStmts is in texture order [checked, yes]
cfg(Block, Entry, Exit) :-
	isBlock(Block), !,
	findall(ChildStmt, hasChild_cfg(Block, ChildStmt), ChildStmts), 
	(% writeln(Block),
		[] = ChildStmts -> Entry = Block, Exit = Block; % empty body
		[First|Rest] = ChildStmts, % otherwise
		cfg(First, FirstEntry, FirstExit),
		% writeln(First),writeln(FirstEntry),
		Entry = FirstEntry, 
		cfg_sq(Rest, FirstExit, Exit)
	).

cfg(Block, Entry, Exit, ParEntry,ParExit,Function) :-
	isBlock(Block), ParEntry == none, !,
	findall(ChildStmt, hasChild_cfg(Block, ChildStmt), ChildStmts), 
	(% writeln(Block),
		[] = ChildStmts -> Entry = Block, Exit = Block; % empty body
		[First|Rest] = ChildStmts, % otherwise
		cfg(First, FirstEntry, FirstExit),
		setBegin(First,Function),
		% writeln(First),writeln(FirstEntry),
		Entry = FirstEntry, 
		cfg_sq(Rest, FirstExit, Exit),
		setEnd(Exit,Function)
	).


cfg(Block, Entry, Exit, ParEntry,ParExit,Function) :-
	isBlock(Block), ParEntry \== none, !,
	findall(ChildStmt, hasChild_cfg(Block, ChildStmt), ChildStmts), 
	(% writeln(Block),
		[] = ChildStmts -> Entry = Block, Exit = Block; % empty body
		[First|Rest] = ChildStmts, % otherwise
		cfg(First, FirstEntry, FirstExit),
		setBegin(ParEntry,Function),
		setNext(ParExit,FirstEntry),
		Entry = FirstEntry, 
		cfg_sq(Rest, FirstExit, Exit),
		setEnd(Exit,Function)
	).



%% The final guard catching the Unknown situation
%% @comment must put at last
cfg(Stmt, Entry, Exit) :-
	c_is_a(Stmt, Class), !,
	write('# Unknown statement: '), write(Stmt), tab(4), writeln(Class),
	Entry = Stmt, Exit = Stmt.



%% helper, recursively process statement sequence
cfg_sq([], Exit, Exit). % last one exit is the exit 
cfg_sq([Stmt|RestStmts], PredExit, Exit) :-
	cfg(Stmt, CurEntry, CurExit),
% write(Stmt), write(' '), write(CurEntry), write(' '), writeln(CurExit),writeln(''),
	(
		PredExit == none -> true; % if previous is break, continue, return, which has exit none, skip
 	setNext(PredExit, CurEntry)
	),
	cfg_sq(RestStmts, CurExit, Exit).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
%
%	Helpers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setBegin(St,Function) :-
	strip_exit(St, P),
	strip_exit(Function, P2),
	writeln(''),writeln(''),
	write(P), write(' pa:beginExecutionInFunction '), write(P2), write(' .'),  nl,
	write(P), write(' pa:firstStatementInFunction '), write(P2), write(' .'),  nl.


setEnd(St,Function) :-
	strip_exit(St, P),
	strip_exit(Function, P2),
	writeln(''),writeln(''),
	write(P), write(' pa:endExecutionInFunction '), write(P2), write(' .'),  nl.



setNext(Pred, Next) :-
	assert(next(Pred, Next)),
%% write(Pred),writeln(Next),
	(
		Pred == none -> true; % if previous is break, continue, return, which has exit none, skip
	dbg_out(Pred, Next, _)
	).

setLastStatement(Last,While) :-
%%	assert(next(Pred, Next, true)),
	dbg_out(Last,While, 'Last').


setTrueNext(Pred, Next) :-
	assert(next(Pred, Next, true)),
	dbg_out(Pred, Next, 'True').

setFalseNext(Pred, Next) :-
	assert(next(Pred, Next, false)),
	dbg_out(Pred, Next, 'False').


%%% get parent structure for break, continue, return
%	@tbd test hasAncestor search behavior

%% get the innermost loop surrounding the break/continue [tested]
getEnclosingLoop(BreakContinueStmt, Loop) :-
	(isBreakStatement(BreakContinueStmt); isContinueStatement(BreakContinueStmt)), !,
	hasAncestor(BreakContinueStmt, Loop), 
	isLoopStatement(Loop), !.

%% get the switch statement of the break [tested]
getEnclosingLoopOrSwitch(BreakStmt, LoopOrSwitchStmt) :-
	isBreakStatement(BreakStmt), !,
	hasAncestor(BreakStmt, LoopOrSwitchStmt),
	(isSwitchStatement(LoopOrSwitchStmt); isLoopStatement(LoopOrSwitchStmt)), !. 

%% get the enclosing function, return function definition
getEnclosingFunction(ReturnStmt, FunctionDef) :-
	isReturnStatement(ReturnStmt), !,
	hasAncestor(ReturnStmt, FunctionDef),
	isFunctionDef(FunctionDef), !.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
%
%	Helpers for debug
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%% node has two type: URL, exit(URL)
strip_exit(ExitNodeInURL, ExitNodeInID) :-
	(
	atom(ExitNodeInURL) -> get_id(ExitNodeInURL, ExitNodeInID_);
		% exit(StmtURL) = ExitNodeInURL,
		% Remove _exit from ExitNodeInURL and bind the remaining to StmtURL
		atomic_list_concat([StmtURL, '_exit'], ExitNodeInURL),
		get_id(StmtURL, Stmt),
		% ExitNodeInID_ = exit(Stmt) dileep
		atomic_list_concat([Stmt, '_exit'], ExitNodeInID_)
	),
%	term_to_atom(ExitNodeInID_, ExitNodeAtom),
	(	
	string_concat(S2,'_exit',ExitNodeInID_) -> write('file0:'), write(ExitNodeInID_), write(' c:hasParent '),write('file0:'), write(S2),write(' .'),nl;
	write('')
	),
	atomic_list_concat(['file0:', ExitNodeInID_], ExitNodeInID).


%% turn off dbg out
%% dbg_out(_, _, _) :- !. 

%  dbg_out(Pred, Next, Label) :-
%	strip_exit(Pred, P), !,
%	strip_exit(Next, N), !,
%	write(P), write(' -> '), write(N), 
%	(
%		nonvar(Label) -> write('[label="'), write(Label), write('"];');
%		write(';')
%	), nl.

dbg_out(Pred, Next, Label) :-
	strip_exit(Pred, P), !,
	strip_exit(Next, N), !,
	
	(nonvar(Label), Label == 'True' -> L2 = ' pa:nextTrueStatement '; 
		(nonvar(Label), Label == 'False' -> L2 = ' pa:nextFalseStatement ';
			(nonvar(Label), Label == 'Last' -> L2 = ' pa:lastStatementInLoop '; 
				(nonvar(Label), Label == 'begin' -> L2 = 'pa:begin'; L2 = ' pa:nextStatement ' 

	)))),

	
	% writeln(L2),
	write(P), write(L2), write(N), write(' .'),
	%( var(Label) -> write(P), write(' -> '), write(N), write(';') ; write('') ),
	%( nonvar(Label), Label == 'True' -> write('[label="'), write(Label), write('"];'), write(';') ; write('') ),
	% ( nonvar(Label), Label == 'False' -> write('[label="'), write(Label), write('"];'), write(';')),
	 nl.


%% end helper for debug


