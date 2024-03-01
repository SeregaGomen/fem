use std::rc::Rc;
use super::error::FemError;

#[derive(Copy, Clone, PartialEq)]
enum Token {
    Number, 
    Sin, 
    Cos, 
    Tan, 
    Exp, 
    Asin, 
    Acos, 
    Atan, 
    Sinh, 
    Cosh, 
    Tanh, 
    Sqrt, 
    Abs, 
    Plus, 
    Minus, 
    Div, 
    Mul, 
    Pow, 
    Or,
    And,
    Not,
    Eq,
    Ne,
    Le,
    Lt,
    Ge,
    Gt,
    Lb,
    Rb,
}

#[derive(Copy, Clone, PartialEq)]
enum TokenType { 
    Delimiter, 
    Numeric, 
    Variable,
    Function, 
    Finished 
}

struct Node {
    val: Rc<f64>,
    tok: Option<Token>,
    children: Vec<Node>,
}

impl Node {
    fn new() -> Self {
        Self{ val: Rc::new(0.0), tok: None, children: vec![] }
    }
    fn double(v: f64) -> Self {
        Self{ val: Rc::new(v), tok: Some(Token::Number), children: vec![] }
    }
    fn unary(t: Token, r: Node) -> Self {
        Self{ val: Rc::new(0.0), tok: Some(t), children: vec![r] }
    }
    fn binary(l: Node, t: Token, r: Node) -> Self {
        Self{ val: Rc::new(0.0), tok: Some(t), children: vec![l, r] }
    }
    fn value(&mut self) -> Result<f64, FemError> {
        match self.tok {
            Some(Token::Number) => Ok(*self.val),
            Some(Token::Plus) => if self.children.len() == 1 { Ok(self.children[0].value()?) } else { Ok(self.children[0].value()? + self.children[1].value()?) }
            Some(Token::Minus) => if self.children.len() == 1 { Ok(-self.children[0].value()?) } else { Ok(self.children[0].value()? - self.children[1].value()?) }
            Some(Token::Mul) => Ok(self.children[0].value()? * self.children[1].value()?),
            Some(Token::Div) => Ok(self.children[0].value()? / self.children[1].value()?),
            Some(Token::Pow) => Ok(self.children[0].value()?.powf(self.children[1].value()?)),
            Some(Token::Sin) => Ok(self.children[0].value()?.sin()),
            Some(Token::Cos) => Ok(self.children[0].value()?.cos()),
            Some(Token::Tan) => Ok(self.children[0].value()?.tan()),
            Some(Token::Abs) => Ok(self.children[0].value()?.abs()),
            Some(Token::Exp) => Ok(self.children[0].value()?.exp()),
            Some(Token::Asin) => Ok(self.children[0].value()?.asin()),
            Some(Token::Acos) => Ok(self.children[0].value()?.acos()),
            Some(Token::Atan) => Ok(self.children[0].value()?.atan()),
            Some(Token::Sinh) => Ok(self.children[0].value()?.sinh()),
            Some(Token::Cosh) => Ok(self.children[0].value()?.cosh()),
            Some(Token::Tanh) => Ok(self.children[0].value()?.tanh()),
            Some(Token::Sqrt) => Ok(self.children[0].value()?.sqrt()),
            Some(Token::And) => if self.children[0].value()? == 1.0 && self.children[1].value()? == 1.0 { Ok(1.0) } else { Ok(0.0) }
            Some(Token::Or) =>  if self.children[0].value()? == 1.0 || self.children[1].value()? == 1.0 { Ok(1.0) } else { Ok(0.0) }
            Some(Token::Not) => if self.children[0].value()? == 1.0 { Ok(0.0) } else { Ok(1.0) }
            Some(Token::Eq) => if self.children[0].value()? == self.children[1].value()? { Ok(1.0) } else { Ok(0.0) }
            Some(Token::Ne) => if self.children[0].value()? != self.children[1].value()? { Ok(1.0) } else { Ok(0.0) }
            Some(Token::Le) => if self.children[0].value()? <= self.children[1].value()? { Ok(1.0) } else { Ok(0.0) }
            Some(Token::Lt) => if self.children[0].value()? < self.children[1].value()? { Ok(1.0) } else { Ok(0.0) }
            Some(Token::Ge) => if self.children[0].value()? >= self.children[1].value()? { Ok(1.0) } else { Ok(0.0) }
            Some(Token::Gt) => if self.children[0].value()? > self.children[1].value()? { Ok(1.0) } else { Ok(0.0) }
            _ => Err(FemError::InternalError),
        }
    }
}

//#[derive(Clone, Debug)]
pub struct Parser<'a> {
    index: usize,
    result: Node,
    variables: Vec<(&'a str, f64)>,
    expression: Vec<char>,
    token: String, 
    tok: Option<Token>,
    token_type: Option<TokenType>,
}

impl<'a> Parser<'a> {
    pub fn new() -> Self {
        Self {
            index: 0,
            result: Node::new(), 
            variables: Vec::new(), 
            expression: Vec::new(), 
            token: String::new(),
            tok: None,
            token_type: None, 
        }
    }
    pub fn set_expression(&mut self, exp: &str) -> Result<(), FemError> {
        self.expression = exp.chars().collect();
        self.compile()
    }
    pub fn set_variable(&mut self, name: &'a str, value: f64) {
        for i in &mut self.variables {
            if i.0 == name {
                i.1 = value;
                return;
            }
        }
        self.variables.push((name, value));
    }
    fn compile(&mut self) -> Result<(), FemError> {
        self.index = 0;
        self.tok = None;
        self.token_type = None;
        loop {
            if self.token_type == Some(TokenType::Finished) {
                break
            }
            self.result = self.get_exp()?;
            if self.token_type == Some(TokenType::Delimiter) {
                if self.tok == Some(Token::Rb) {
                    return Err(FemError::BracketError);
                } else {
                    return Err(FemError::SyntaxError);
                }
            }
        }
        Ok(())
    }
    pub fn value(&mut self) -> Result<f64, FemError> {
        self.result.value()
    }
    fn is_expression(&self) -> bool {
        if self.index < self.expression.len() {
            return true
        }   
        false 
    }
    fn get_token(&mut self) -> Result<Option<TokenType>, FemError> {
        self.token = String::new();
        self.token_type = None;
        self.tok = None;
        // Обработка пустой строки
        if !self.is_expression() {
            self.token_type = Some(TokenType::Finished);
            return Ok(self.token_type)
        }
        // Пропуск ведущих пробелов
        for i in self.index..self.expression.len() {
            if self.expression[i] != ' ' && self.expression[i] != '\t' {
                break
            }
            self.index += 1
        }
        // Обработка разделителя
        if self.is_expression() && "+-*/()=^<>!".contains(self.expression[self.index]) {
            self.token.push(self.expression[self.index]);
            self.index += 1;
            // Проверка на наличие двойного разделителя
            if self.is_expression() && "=<>".contains(self.expression[self.index]) {
                self.token.push(self.expression[self.index]);
                self.index += 1;
            }
            if !self.find_delimiter() {
                return Err(FemError::SyntaxError);       
            }
            self.token_type = Some(TokenType::Delimiter);
            return Ok(self.token_type)
        }
        // Обработка числа
        if self.is_expression() && self.expression[self.index].is_ascii_digit() {
            while self.is_expression() && self.expression[self.index].is_ascii_digit() {
                self.token.push(self.expression[self.index]);
                self.index += 1
            }
            if self.is_expression() && self.expression[self.index] == '.' {
                self.token.push('.');
                self.index += 1;
                while self.is_expression() && self.expression[self.index].is_ascii_digit() {
                    self.token.push(self.expression[self.index]);
                    self.index += 1
                }
            }
            if self.is_expression() && (self.expression[self.index].to_ascii_uppercase() == 'E') {
                self.token.push('E');
                self.index += 1;
                if self.is_expression() && (self.expression[self.index] == '+' || self.expression[self.index] == '-') {
                    self.token.push(self.expression[self.index]);
                    self.index += 1;
                    while self.is_expression() && self.expression[self.index].is_ascii_digit() {
                        self.token.push(self.expression[self.index]);
                        self.index += 1
                    }
                } else {
                    return Err(FemError::InvalidNumber);   
                }
            }
            self.token_type = Some(TokenType::Numeric);
            return Ok(self.token_type)
        }
        // Обработка строкового литерала
        if self.index < self.expression.len() && (self.expression[self.index].is_ascii_alphabetic() || self.expression[self.index] == '_') {
            while self.index < self.expression.len() && (self.expression[self.index].is_ascii_alphabetic() || self.expression[self.index].is_ascii_digit() || self.expression[self.index] == '_') {
                self.token.push(self.expression[self.index]);
                self.index += 1
            }   
            self.token_type = None;
            if !self.find_delimiter() {
                if !self.find_variable() {
                    if !self.find_function() {
                        return Err(FemError::UndefError);
                    }
                }
            }   
            return Ok(self.token_type);
        }
        Err(FemError::SyntaxError)
    }
    fn find_delimiter(&mut self) -> bool {
        self.tok = match &self.token[..] {
            "+" => Some(Token::Plus),
            "-" => Some(Token::Minus),
            "/" => Some(Token::Div),
            "*" => Some(Token::Mul),
            "^" => Some(Token::Pow),
            "or" => Some(Token::Or),
            "and" => Some(Token::And),
            "not" => Some(Token::Not),
            "==" => Some(Token::Eq),
            "!=" => Some(Token::Ne),
            "<=" => Some(Token::Le),
            "<" => Some(Token::Lt),
            ">=" => Some(Token::Ge),
            ">" => Some(Token::Gt),
            "(" => Some(Token::Lb),
            ")" => Some(Token::Rb),
            _ => return false,
        };
        self.token_type = Some(TokenType::Delimiter);
        true
    }
    fn find_function(&mut self) -> bool {
        self.tok = match &self.token[..] {
            "sqrt" => Some(Token::Sqrt),
            "sin" => Some(Token::Sin),
            "cos" => Some(Token::Cos),
            "tan" => Some(Token::Tan),
            "exp" => Some(Token::Exp),
            "asin" => Some(Token::Asin),
            "acos" => Some(Token::Acos),
            "atan" => Some(Token::Atan),
            "sinh" => Some(Token::Sinh),
            "cosh" => Some(Token::Cosh),
            "tanh" => Some(Token::Tanh),
            "abs" => Some(Token::Abs),
            _ => return false,
        };
        self.token_type = Some(TokenType::Function);
        true
    }
    fn find_variable(&mut self) ->bool {
        let mut res = false;
        for i in &self.variables {
            if i.0 == self.token {
                self.token_type = Some(TokenType::Variable);
                res = true;
                break;
            }
        }
        res
    }
    fn get_exp(&mut self) -> Result<Node, FemError> {
        self.get_token()?;
        self.token_or()
    }
    fn token_or(&mut self) -> Result<Node, FemError> {
        let mut res = self.token_and()?;
        while self.token_type != Some(TokenType::Finished) && self.tok == Some(Token::Or) {
            self.get_token()?;
            let hold = self.token_and()?;
            res = Node::binary(res, Token::Or, hold);
        }
        Ok(res)   
    }
    fn token_and(&mut self) -> Result<Node, FemError> {
        let mut res = self.token_not()?;
        while self.token_type != Some(TokenType::Finished) && self.tok == Some(Token::And) {
            self.get_token()?;
            let hold = self.token_not()?;
            res = Node::binary(res, Token::And, hold);
        }
        Ok(res)   
    }
    fn token_not(&mut self) -> Result<Node, FemError> {
        let op: Option<Token> = if self.tok == Some(Token::Not) { Some(Token::Not) } else { None }; 
        if self.token_type != Some(TokenType::Finished) && self.tok == Some(Token::Not) {
            self.get_token()?;    
        }
        let mut res = self.token_eq()?;
        if op == Some(Token::Not) {
            res = Node::unary(Token::Not, res);
        } 
        Ok(res)   
    }
    fn token_eq(&mut self) -> Result<Node, FemError> {
        let mut res = self.token_add()?;
        while self.token_type != Some(TokenType::Finished) && (self.tok == Some(Token::Gt) || self.tok == Some(Token::Ge) || 
            self.tok == Some(Token::Lt) || self.tok == Some(Token::Le) || self.tok == Some(Token::Eq) || self.tok == Some(Token::Ne)) {
            let op = self.tok.unwrap(); 
            self.get_token()?;
            let hold = self.token_add()?;
            res = Node::binary(res, op, hold);
        }
        Ok(res)   
    }
    fn token_add(&mut self) -> Result<Node, FemError> {
        let mut res = self.token_mul()?;
        while self.token_type != Some(TokenType::Finished) && (self.tok == Some(Token::Plus) || self.tok == Some(Token::Minus)) {
            let op = self.tok.unwrap(); 
            self.get_token()?;
            let hold = self.token_mul()?;
            res = Node::binary(res, op, hold);
        }
        Ok(res)   
    }
    fn token_mul(&mut self) -> Result<Node, FemError> {
        let mut res = self.token_pow()?;
        while self.token_type != Some(TokenType::Finished) && (self.tok == Some(Token::Mul) || self.tok == Some(Token::Div)) {
            let op: Token = if self.tok == Some(Token::Mul) { Token::Mul } else { Token::Div }; 
            self.get_token()?;
            let hold = self.token_pow()?;
            res = Node::binary(res, op, hold);
        }
        Ok(res)   
    }
    fn token_pow(&mut self) -> Result<Node, FemError> {
        let mut res = self.token_un()?;
        if self.token_type != Some(TokenType::Finished) && self.tok == Some(Token::Pow) {
            self.get_token()?;
            let hold = self.token_bracket()?;
            res = Node::binary(res, Token::Pow, hold);
        }
        Ok(res)   
    }
    fn token_un(&mut self) -> Result<Node, FemError> {
        let op: Option<Token> = if self.tok == Some(Token::Plus) { Some(Token::Plus) } else { if self.tok == Some(Token::Minus) { Some(Token::Minus) } else { None } };
        if self.token_type == Some(TokenType::Delimiter) && (self.tok == Some(Token::Plus) || self.tok == Some(Token::Minus)) {
            self.get_token()?;
        }
        let mut res = self.token_bracket()?;
        if op != None {
            res = Node::unary(op.unwrap(), res);
        }
        Ok(res)
    }
    fn token_bracket(&mut self) -> Result<Node, FemError> {
        let res;
        if self.tok == Some(Token::Lb) && self.token_type == Some(TokenType::Delimiter) {
            self.get_token()?;
            res = self.token_or()?;
            if self.tok != Some(Token::Rb) {
                return Err(FemError::SyntaxError);
            }
            self.get_token()?;
        } else {
            res = self.token_prim()?;
        }
        Ok(res)
    }
    fn token_prim(&mut self) -> Result<Node, FemError> {
        let mut is_find = false;
        let mut res = Err(FemError::SyntaxError);
        match self.token_type {
            Some(TokenType::Numeric) => {
                res = Ok(Node::double(self.token.parse().unwrap()));
                self.get_token()?;
            }
            Some(TokenType::Variable) => {
                for i in 0..self.variables.len() {
                    if self.variables[i].0 == self.token {
                        res = Ok(Node::double(self.variables[i].1));
                        self.get_token()?;
                        is_find = true;
                        break;
                    }
                }
                if !is_find {
                    res = Err(FemError::UndefError);
                }
            }
            Some(TokenType::Function) => res = Ok(self.token_func()?),
            _ => res = Err(FemError::SyntaxError),
        }
        //self.get_token()?;
        res
    }
    fn token_func(&mut self) -> Result<Node, FemError> {
        let mut res;
        let fun_tok = self.tok.unwrap();
        self.get_token()?;
        if self.token.len() == 0 || self.tok != Some(Token::Lb) {
            return Err(FemError::SyntaxError);
        }
        self.get_token()?;
        res = self.token_add()?;
        res = Node::unary(fun_tok, res);
        if self.tok != Some(Token::Rb) {
            return Err(FemError::SyntaxError);     
        }
        self.get_token()?;
        Ok(res)
    }
}
