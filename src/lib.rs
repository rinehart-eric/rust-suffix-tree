extern crate bio;

use std::u8;
use std::i32;
use std::rc::Rc;
use std::cell::RefCell;
use std::collections::HashMap;
use std::borrow::Borrow;
use std::borrow::BorrowMut;
use self::bio::alphabets::Alphabet;

/// A struct representing a transition (or edge) from one explicit node to another
#[derive(Debug)]
struct Transition {
    /// The edge's string slice (inclusive)
    slice: (i32, i32),
    /// The index of the further node from the root
    next: usize,
}

/// An explicit node in the suffix tree
#[derive(Debug)]
struct Node {
    /// `Transition`s mapped by the first character along each edge
    transitions: HashMap<u8, Transition>,
    /// Optional suffix link
    suffix_link: Option<usize>,
}

impl Node {
    /// Creates a new explicit `Node` with an empty transition map and no suffix link
    pub fn new() -> Node {
        Node {
            transitions: HashMap::new(),
            suffix_link: None,
        }
    }
}

/// A struct which implements Ukkonen's suffix tree construction algorithm.
///
/// # Examples
///
/// ```
/// let str = "ATCG".to_string();
/// let tree = suffix_tree::SuffixTree::new(bio::alphabets::dna::alphabet(), &str);
/// let search = "TC".to_string();
/// if tree.contains(&search) {
///     println!("{} contains {}!", str, search);
/// }
/// ```
#[derive(Debug)]
pub struct SuffixTree {
    nodes: Vec<Rc<RefCell<Node>>>,
    string: Vec<u8>,
}

impl SuffixTree {
    /// Creates a new suffix tree from `string` constrained by `alphabet`.
    pub fn new(alphabet: Alphabet, string: &String) -> SuffixTree {
        let mut tree_string = string.chars().map(|c| c as u8).collect::<Vec<u8>>();
        if !alphabet.is_word(&tree_string) {
            panic!("String is not constrained to the given alphabet!");
        }
        let end_char = '$' as u8;
        tree_string.push(end_char);
        let mut new_tree = SuffixTree {
            nodes: Vec::new(),
            string: tree_string,
        };

        // Node 0: auxiliary node (allows root node to not be a special case)
        let auxiliary = Rc::new(RefCell::new(Node::new()));
        new_tree.nodes.push(auxiliary.clone());

        // Node 1: root node
        let root = Rc::new(RefCell::new(Node::new()));
        new_tree.nodes.push(root.clone());

        for (i, sym) in alphabet.symbols.into_iter().enumerate() {
            let t = sym as u8;
            let idx = -(i as i32) - 2;
            (*auxiliary).borrow_mut().transitions.insert(t, Transition {
                slice: (idx, idx),
                next: 1,
            });
        }
        (*auxiliary).borrow_mut().transitions.insert(end_char, Transition {
            slice: (-1, -1),
            next: 1,
        });
        (*root).borrow_mut().suffix_link = Some(0);

        new_tree.build();
        new_tree
    }

    fn build(&mut self) {
        let mut s = 1; // root node
        let mut k = 1;
        let mut i = 0;
        for _ in &self.string {
            i += 1;
            match SuffixTree::update(&mut self.nodes, &self.string, s, (k, i as i32)) {
                (s_p, k_p) => { s = s_p; k = k_p; }
            }
            match SuffixTree::canonize(&self.nodes, &self.string, s, (k, i as i32)) {
                (s_p, k_p) => { s = s_p; k = k_p; }
            }
        }

        let len = self.string.len() as i32;
        for node_ref in &mut self.nodes {
            let mut node = (**node_ref).borrow_mut();
            for (_, trans) in &mut node.transitions {
                if trans.slice.1 == i32::MAX {
                    trans.slice.1 = len;
                }
            }
        }
    }

    fn update(nodes: &mut Vec<Rc<RefCell<Node>>>, string: &Vec<u8>, node: usize, pair: (i32, i32)) -> (usize, i32) {
        let (mut k, i) = pair;
        let mut s = node;

        let mut oldr = 1;

        let (mut endpoint, mut r) = SuffixTree::test_and_split(nodes, string, s, (k, i - 1), SuffixTree::char_at(string, i));
        while !endpoint {
            let t_i = SuffixTree::char_at(string, i);
            let r_p = Rc::new(RefCell::new(Node::new()));
            let r_p_id = nodes.len();
            nodes.push(r_p);

            let trans = Transition {
                slice: (i, i32::MAX),
                next: r_p_id,
            };

            let rc = nodes[r].clone();
            (*rc).borrow_mut().transitions.insert(t_i, trans);

            if oldr > 1 { // not root or aux
                let oldr_refcell = nodes[oldr].clone();
                (*oldr_refcell).borrow_mut().suffix_link = Some(r);
            }
            oldr = r;
            {
                let s_ref = nodes[s].clone();
                let s_node = (*s_ref).borrow();
                match s_node.suffix_link.map(|link| SuffixTree::canonize(nodes, string, link, (k, i - 1))) {
                    Some((s_p, k_p)) => {
                        s = s_p;
                        k = k_p;
                    },
                    None => panic!("No suffix link for {:?}\nin {:?}", s_node, nodes),
                }
            }

            let (new_endpoint, new_r) = SuffixTree::test_and_split(nodes, string, s, (k, i - 1), SuffixTree::char_at(string, i));
            endpoint = new_endpoint;
            r = new_r;
        }
        if oldr > 1 { // not root or aux
            let oldr_refcell = nodes[oldr].clone();
            (*oldr_refcell).borrow_mut().suffix_link = Some(s);
        }
        (s, k)
    }

    fn canonize(nodes: &Vec<Rc<RefCell<Node>>>, string: &Vec<u8>, s: usize, pair: (i32, i32)) -> (usize, i32) {
        let (mut k, p) = pair;
        if p < k {
            return (s, k);
        }

        let ((mut k_p, mut p_p), mut s_p): ((i32, i32), usize);

        let mut curr = (nodes[s].clone(), s);
        let mut t_k = SuffixTree::char_at(string, k);
        match (*curr.0).borrow().transitions.get(&t_k) {
            Some(t) => {
                k_p = t.slice.0;
                p_p = t.slice.1;
                s_p = t.next;
            },
            None => panic!("No {}-transition for {:?}\nin {:?}", t_k, (*curr.0).borrow(), nodes),
        };
        while (p_p - k_p) <= (p - k) {
            k += p_p - k_p + 1;
            curr = (nodes[s_p].clone(), s_p);

            if k <= p {
                let curr_node = (*curr.0).borrow();
                t_k = SuffixTree::char_at(string, k);
                match curr_node.transitions.get(&t_k) {
                    Some(t) => {
                        k_p = t.slice.0;
                        p_p = t.slice.1;
                        s_p = t.next;
                    },
                    None => panic!("No {}-transition for {:?}\nin {:?}", t_k, curr_node, nodes),
                };
            }
        }
        (curr.1, k)
    }

    fn test_and_split(nodes: &mut Vec<Rc<RefCell<Node>>>, string: &Vec<u8>, s: usize, pair: (i32, i32), t: u8) -> (bool, usize) {
        let (k, p) = pair;
        let s_ref = nodes[s].clone();
        let mut s_node = (*s_ref).borrow_mut();

        if k > p {
            return (s_node.borrow().transitions.contains_key(&t), s);
        }

        let trans_char = SuffixTree::char_at(string, k);
        let trans = match s_node.transitions.get_mut(&trans_char) {
            Some(t) => t,
            _ => panic!("No {}-transition for {:?}", trans_char, s),
        };
        let (k_p, p_p) = trans.slice;

        if t == SuffixTree::char_at(string, (k_p + p - k + 1)) {
            return (true, s);
        } else {
            // Create the new split node
            let r = Rc::new(RefCell::new(Node::new()));
            let r_id = nodes.len();
            nodes.push(r);

            let r_k = k_p + p - k + 1;
            let r_to_s_p = Transition {
                slice: (r_k, p_p),
                next: trans.next
            };

            // Turn the s -> s' transition into the s -> r transition
            *trans = Transition {
                slice: (k_p, (k_p + p - k)),
                next: r_id,
            };

            // Insert the r -> s' transition into r's transition map
            let r_ref = nodes[r_id].clone();
            let mut r_node = (*r_ref).borrow_mut();
            r_node.transitions.insert(SuffixTree::char_at(string, r_k), r_to_s_p);

            return (false, r_id);
        }
    }

    /// Returns the `u8` at position (`idx` - 1) in the string (i.e., 1-indexed position)
    fn char_at(str: &Vec<u8>, idx: i32) -> u8 {
        str[(idx - 1) as usize]
    }

    /// Returns `true` if this suffix tree contains `str`, and `false` if not.
    pub fn contains(&self, str: &String) -> bool {
        let mut curr = 1;
        let mut chars = str.chars();

        loop {
            let c = chars.next().map(|c| c as u8);
            if c.is_none() {
                return true;
            }
            let curr_node = (*self.nodes[curr]).borrow();
            match curr_node.transitions.get(&c.unwrap()) {
                Some(trans) => {
                    for i in trans.slice.0..trans.slice.1 {
                        match chars.next().map(|ch| ch as u8) {
                            Some(ch) => {
                                if ch != SuffixTree::char_at(&self.string, i + 1) {
                                    return false;
                                }
                            },
                            None => return true,
                        }
                    }
                    curr = trans.next;
                },
                None => return false,
            }
        }
    }
}
