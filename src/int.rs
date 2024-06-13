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

impl Character for u32 {
    #[inline]
    fn char_to_usize(self) -> usize {
        self as usize
    }

    #[inline]
    fn zero() -> u32 {
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
    #[inline(always)]
    fn as_usize(self) -> usize {
        self.to_usize().unwrap()
    }
}

impl Index for i32 { }
impl Index for i64 { }
impl Index for i16 { }
// impl Index for i8  { }

pub trait Counter: Index {
    fn is_empty(&self) -> bool;

    fn is_counter(&self) -> bool;

    fn get_count(&self) -> Self;

    fn inc_counter(&mut self) {
        *self -= Self::ONE;
    }

    fn init_counter(&mut self) {
        *self = Self::max_value();
    }
}
