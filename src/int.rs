use num_traits::*;

pub trait Character: Copy + Ord + std::fmt::Debug {
    fn char_to_usize(self) -> usize;

    fn zero() -> Self;
}

impl Character for i32 {
    #[inline]
    fn char_to_usize(self) -> usize {
        self as usize
    }

    #[inline]
    fn zero() -> i32 {
        0
    }
}

impl Character for i64 {
    #[inline]
    fn char_to_usize(self) -> usize {
        self as usize
    }

    #[inline]
    fn zero() -> i64 {
        0
    }
}

impl Character for i16 {
    #[inline]
    fn char_to_usize(self) -> usize {
        self as usize
    }

    #[inline]
    fn zero() -> i16 {
        0
    }

}

impl Character for u8 {
    #[inline]
    fn char_to_usize(self) -> usize {
        self as usize
    }

    #[inline]
    fn zero() -> u8 {
        0
    }
}


pub trait Index: Signed + PrimInt + NumAssign + ConstZero + ConstOne + NumCast + Bounded + Character + Copy {
    fn as_usize(self) -> usize {
        self.to_usize().unwrap()
    }
}

impl Index for i32 { }
impl Index for i64 { }
impl Index for i16 { }
// impl Index for i8  { }
