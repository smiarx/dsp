	.file	"Reverb.cpp"
	.text
	.p2align 4
	.globl	_Z13SCReverb_DtorP8SCReverb
	.type	_Z13SCReverb_DtorP8SCReverb, @function
_Z13SCReverb_DtorP8SCReverb:
.LFB2338:
	.cfi_startproc
	movq	_ZL2ft(%rip), %rax
	pushq	%rbx
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	movq	%rdi, %rbx
	movq	200(%rdi), %rsi
	movq	(%rdi), %rdi
	call	*128(%rax)
	movq	_ZL2ft(%rip), %rax
	movq	176(%rbx), %rsi
	movq	(%rbx), %rdi
	call	*128(%rax)
	movq	_ZL2ft(%rip), %rax
	movq	192(%rbx), %rsi
	movq	(%rbx), %rdi
	popq	%rbx
	.cfi_def_cfa_offset 8
	movq	128(%rax), %rax
	jmp	*%rax
	.cfi_endproc
.LFE2338:
	.size	_Z13SCReverb_DtorP8SCReverb, .-_Z13SCReverb_DtorP8SCReverb
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC3:
	.string	"SCReverb_Ctor"
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC4:
	.string	"%s: alloc failed, increase server's RT memory (e.g. via ServerOptions)\n"
	.text
	.p2align 4
	.globl	_Z13SCReverb_CtorP8SCReverb
	.type	_Z13SCReverb_CtorP8SCReverb, @function
_Z13SCReverb_CtorP8SCReverb:
.LFB2328:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	leaq	_Z13SCReverb_nextP8SCReverbi(%rip), %rax
	pxor	%xmm0, %xmm0
	pxor	%xmm1, %xmm1
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	movl	$1024, %esi
	movq	%rdi, %rbx
	subq	$120, %rsp
	.cfi_def_cfa_offset 144
	movaps	.LC0(%rip), %xmm2
	movq	%fs:40, %rbp
	movq	%rbp, 104(%rsp)
	xorl	%ebp, %ebp
	movq	%rax, 88(%rdi)
	movq	_ZL2ft(%rip), %rax
	movaps	%xmm2, (%rsp)
	movaps	.LC1(%rip), %xmm2
	movq	$0, 200(%rdi)
	movaps	%xmm2, 16(%rsp)
	movaps	.LC2(%rip), %xmm2
	movups	%xmm1, 152(%rdi)
	movaps	%xmm2, 32(%rsp)
	movdqa	.LC8(%rip), %xmm2
	movups	%xmm0, 168(%rdi)
	movups	%xmm2, 104(%rdi)
	movdqa	.LC9(%rip), %xmm2
	movups	%xmm0, 184(%rdi)
	movups	%xmm2, 120(%rdi)
	movdqa	.LC10(%rip), %xmm2
	movq	$0, 96(%rsp)
	movups	%xmm2, 136(%rdi)
	movq	(%rdi), %rdi
	movaps	%xmm0, 64(%rsp)
	movaps	%xmm0, 80(%rsp)
	movaps	%xmm1, 48(%rsp)
	call	*112(%rax)
	movq	_ZL2ft(%rip), %rdx
	leaq	.LC3(%rip), %rsi
	leaq	.LC4(%rip), %rdi
	testq	%rax, %rax
	je	.L11
	movq	%rax, 200(%rbx)
	movq	(%rbx), %rdi
	movl	$8192, %esi
	call	*112(%rdx)
	movq	_ZL2ft(%rip), %r8
	movq	%rax, %rdx
	testq	%rax, %rax
	je	.L12
	leaq	8(%rax), %rdi
	movq	%rax, %rcx
	movq	$0, (%rax)
	movl	$262144, %esi
	movq	$0, 5336(%rax)
	andq	$-8, %rdi
	movq	%rbp, %rax
	subq	%rdi, %rcx
	addl	$5344, %ecx
	shrl	$3, %ecx
	rep stosq
	movq	%rdx, 176(%rbx)
	movq	(%rbx), %rdi
	call	*112(%r8)
	testq	%rax, %rax
	je	.L13
	movl	$143048, %edx
	xorl	%esi, %esi
	movq	%rax, %rdi
	call	memset@PLT
	movq	%rax, 192(%rbx)
	movq	80(%rbx), %rax
	movq	(%rax), %rdx
	movq	8(%rax), %rax
	movl	$0x00000000, (%rdx)
	movl	$0x00000000, (%rax)
.L4:
	movq	104(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L14
	addq	$120, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L13:
	.cfi_restore_state
	movq	_ZL2ft(%rip), %rdx
	leaq	.LC3(%rip), %rsi
	leaq	.LC4(%rip), %rdi
	xorl	%eax, %eax
.L11:
	call	*32(%rdx)
	movq	_ZL2ft(%rip), %rax
	movq	80(%rax), %rax
	movq	%rax, 88(%rbx)
	movl	$1, %eax
	movw	%ax, 38(%rbx)
	jmp	.L4
	.p2align 4,,10
	.p2align 3
.L12:
	leaq	.LC3(%rip), %rsi
	leaq	.LC4(%rip), %rdi
	xorl	%eax, %eax
	call	*32(%r8)
	movq	_ZL2ft(%rip), %rax
	movl	$1, %edx
	movq	80(%rax), %rax
	movw	%dx, 38(%rbx)
	movq	%rax, 88(%rbx)
	jmp	.L4
.L14:
	call	__stack_chk_fail@PLT
	.cfi_endproc
.LFE2328:
	.size	_Z13SCReverb_CtorP8SCReverb, .-_Z13SCReverb_CtorP8SCReverb
	.align 2
	.p2align 4
	.globl	_ZN6Reverb7processEPPfS1_i
	.type	_ZN6Reverb7processEPPfS1_i, @function
_ZN6Reverb7processEPPfS1_i:
.LFB2309:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	movq	96(%rdi), %rax
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	movl	80(%rdi), %r15d
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movl	%ecx, %ebp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movq	(%rsi), %r9
	movq	8(%rsi), %r10
	movl	64(%rdi), %esi
	movq	(%rdx), %r13
	movq	8(%rdx), %r14
	movl	%ecx, -104(%rsp)
	movl	$128, %ecx
	movq	88(%rdi), %rdx
	movl	%esi, -116(%rsp)
	movq	72(%rdi), %rsi
	cmpl	%ecx, %ebp
	cmovle	%ebp, %ecx
	movq	%rsi, -16(%rsp)
	testl	%ebp, %ebp
	je	.L16
	movslq	%ecx, %r8
	movl	$4, %r11d
	movl	%ebp, -52(%rsp)
	movq	.LC12(%rip), %xmm2
	leaq	0(,%r8,8), %rbx
	movl	%ecx, %r8d
	movq	.LC13(%rip), %xmm1
	movq	%rbx, -48(%rsp)
	leaq	0(,%r8,4), %rbx
	cmovle	%r11, %rbx
	salq	$3, %r8
	testl	%ebp, %ebp
	movq	%rbx, -112(%rsp)
	leaq	-4(%r8), %rbx
	cmovg	%rbx, %r11
	movq	%r11, -32(%rsp)
	movl	$1, %r11d
	cmovg	%ecx, %r11d
	xorl	%ebx, %ebx
	movl	%r11d, %r12d
	movl	%r11d, -100(%rsp)
	leal	-1(%rcx), %r11d
	salq	$2, %r11
	testl	%ebp, %ebp
	cmovg	%r11, %rbx
	movl	$8, %r11d
	cmovg	%r8, %r11
	movl	%r12d, %r8d
	andl	$-4, %r12d
	movq	%rbx, -40(%rsp)
	movl	%r12d, %ebx
	shrl	$2, %r8d
	movq	%rbx, -88(%rsp)
	salq	$3, %rbx
	salq	$4, %r8
	movq	%r11, -24(%rsp)
	movl	%r12d, -56(%rsp)
	movq	%rbx, -80(%rsp)
	movl	%r15d, %ebx
	movq	%r14, %r15
	movl	%ebx, %r14d
	.p2align 4
	.p2align 3
.L44:
	movl	-104(%rsp), %ebx
	testl	%ebx, %ebx
	jle	.L18
	cmpl	$3, %ebx
	jle	.L46
	movq	-112(%rsp), %r12
	movq	-32(%rsp), %rbx
	leaq	(%r9,%r12), %r11
	addq	%rax, %rbx
	cmpq	%r11, %rax
	setnb	%bpl
	cmpq	%rbx, %r9
	setnb	%r11b
	orb	%bpl, %r11b
	je	.L46
	leaq	(%r10,%r12), %r11
	cmpq	%r11, %rax
	setnb	%r11b
	cmpq	%rbx, %r10
	setnb	%bl
	orb	%r11b, %bl
	je	.L46
	xorl	%r11d, %r11d
	.p2align 4
	.p2align 3
.L20:
	movups	(%r10,%r11), %xmm0
	movups	(%r9,%r11), %xmm7
	addps	%xmm7, %xmm0
	movaps	%xmm0, %xmm3
	movss	%xmm0, (%rax,%r11,2)
	shufps	$85, %xmm0, %xmm3
	movss	%xmm3, 8(%rax,%r11,2)
	movaps	%xmm0, %xmm3
	unpckhps	%xmm0, %xmm3
	shufps	$255, %xmm0, %xmm0
	movss	%xmm0, 24(%rax,%r11,2)
	movss	%xmm3, 16(%rax,%r11,2)
	addq	$16, %r11
	cmpq	%r8, %r11
	jne	.L20
	movq	-88(%rsp), %r11
	movq	-80(%rsp), %rbx
	salq	$2, %r11
	addq	%rax, %rbx
	leaq	(%r10,%r11), %rbp
	addq	%r9, %r11
	testb	$3, -100(%rsp)
	je	.L23
	movss	(%r11), %xmm0
	movl	-56(%rsp), %r12d
	addss	0(%rbp), %xmm0
	addl	$1, %r12d
	movss	%xmm0, (%rbx)
	cmpl	%r12d, %ecx
	jle	.L23
	movss	4(%rbp), %xmm0
	movl	-56(%rsp), %r12d
	addss	4(%r11), %xmm0
	addl	$2, %r12d
	movss	%xmm0, 8(%rbx)
	cmpl	%r12d, %ecx
	jle	.L23
	movss	8(%rbp), %xmm0
	addss	8(%r11), %xmm0
	movss	%xmm0, 16(%rbx)
.L23:
	movq	-40(%rsp), %rbx
	movl	-116(%rsp), %r11d
	movss	.LC11(%rip), %xmm0
	addq	$4, %rbx
	movq	%rbx, -96(%rsp)
	addq	%rbx, %r10
	addq	%rbx, %r9
	xorl	%ebx, %ebx
	.p2align 4
	.p2align 3
.L24:
	leal	150(%r11), %ebp
	movss	(%rax,%rbx,8), %xmm3
	andl	$2047, %ebp
	movss	(%rsi,%rbp,4), %xmm5
	movl	%r11d, %ebp
	subl	$1, %r11d
	andl	$2047, %ebp
	andl	$2047, %r11d
	movaps	%xmm5, %xmm4
	mulss	%xmm0, %xmm4
	subss	%xmm4, %xmm3
	movaps	%xmm3, %xmm4
	mulss	%xmm0, %xmm4
	addss	%xmm5, %xmm4
	movss	%xmm4, (%rax,%rbx,8)
	addq	$1, %rbx
	movss	%xmm3, (%rsi,%rbp,4)
	cmpl	%ebx, %ecx
	jg	.L24
	movl	-116(%rsp), %r11d
	xorl	%ebx, %ebx
	.p2align 4
	.p2align 3
.L25:
	leal	363(%r11), %ebp
	movss	(%rax,%rbx,8), %xmm3
	andl	$2047, %ebp
	movss	(%rsi,%rbp,4), %xmm5
	leal	150(%r11), %ebp
	subl	$1, %r11d
	andl	$2047, %ebp
	andl	$2047, %r11d
	movaps	%xmm5, %xmm4
	mulss	%xmm0, %xmm4
	subss	%xmm4, %xmm3
	movaps	%xmm3, %xmm4
	mulss	%xmm0, %xmm4
	addss	%xmm5, %xmm4
	movss	%xmm4, (%rax,%rbx,8)
	addq	$1, %rbx
	movss	%xmm3, (%rsi,%rbp,4)
	cmpl	%ebx, %ecx
	jg	.L25
	movl	-116(%rsp), %r11d
	xorl	%ebx, %ebx
	.p2align 4
	.p2align 3
.L26:
	leal	682(%r11), %ebp
	movss	(%rax,%rbx,8), %xmm3
	andl	$2047, %ebp
	movss	(%rsi,%rbp,4), %xmm5
	leal	363(%r11), %ebp
	subl	$1, %r11d
	andl	$2047, %ebp
	andl	$2047, %r11d
	movaps	%xmm5, %xmm4
	mulss	%xmm0, %xmm4
	subss	%xmm4, %xmm3
	movaps	%xmm3, %xmm4
	mulss	%xmm0, %xmm4
	addss	%xmm5, %xmm4
	movss	%xmm4, (%rax,%rbx,8)
	addq	$1, %rbx
	movss	%xmm3, (%rsi,%rbp,4)
	cmpl	%ebx, %ecx
	jg	.L26
	movl	-116(%rsp), %r11d
	xorl	%ebx, %ebx
	.p2align 4
	.p2align 3
.L27:
	leal	1208(%r11), %ebp
	movss	(%rax,%rbx,8), %xmm3
	andl	$2047, %ebp
	movss	(%rsi,%rbp,4), %xmm5
	leal	682(%r11), %ebp
	subl	$1, %r11d
	andl	$2047, %ebp
	andl	$2047, %r11d
	movaps	%xmm5, %xmm4
	mulss	%xmm0, %xmm4
	subss	%xmm4, %xmm3
	movaps	%xmm3, %xmm4
	mulss	%xmm0, %xmm4
	addss	%xmm5, %xmm4
	movss	%xmm4, (%rax,%rbx,8)
	addq	$1, %rbx
	movss	%xmm3, (%rsi,%rbp,4)
	cmpl	%ebx, %ecx
	jg	.L27
	xorl	%r11d, %r11d
	.p2align 5
	.p2align 4
	.p2align 3
.L28:
	movss	(%rax,%r11,8), %xmm0
	movss	%xmm0, 4(%rax,%r11,8)
	addq	$1, %r11
	cmpl	%r11d, %ecx
	jg	.L28
	movl	%r14d, %r11d
	xorl	%ebx, %ebx
	.p2align 4
	.p2align 3
.L29:
	leal	15952(%r11), %ebp
	leal	17753(%r11), %r12d
	subl	$1, %r11d
	andl	$32767, %ebp
	andl	$32767, %r12d
	andl	$32767, %r11d
	movss	(%rdx,%r12,8), %xmm3
	movss	4(%rdx,%rbp,8), %xmm0
	unpcklps	%xmm3, %xmm0
	movq	(%rax,%rbx,8), %xmm3
	mulps	%xmm2, %xmm0
	addps	%xmm3, %xmm0
	movlps	%xmm0, (%rax,%rbx,8)
	addq	$1, %rbx
	cmpl	%ebx, %ecx
	jg	.L29
	movss	52(%rdi), %xmm4
	movss	48(%rdi), %xmm3
	xorl	%r11d, %r11d
	.p2align 4
	.p2align 3
.L30:
	movss	56(%rdi), %xmm6
	movss	40(%rdi), %xmm7
	movaps	%xmm3, %xmm0
	movss	(%rax,%r11,8), %xmm3
	movss	60(%rdi), %xmm5
	mulss	%xmm6, %xmm7
	mulss	16(%rdi), %xmm6
	subss	%xmm7, %xmm3
	movss	32(%rdi), %xmm7
	mulss	%xmm0, %xmm7
	mulss	8(%rdi), %xmm0
	subss	%xmm7, %xmm3
	addss	%xmm6, %xmm0
	movss	(%rdi), %xmm6
	mulss	%xmm3, %xmm6
	addss	%xmm6, %xmm0
	movss	%xmm0, (%rax,%r11,8)
	movss	44(%rdi), %xmm6
	movaps	%xmm4, %xmm0
	movss	4(%rax,%r11,8), %xmm4
	mulss	%xmm5, %xmm6
	mulss	20(%rdi), %xmm5
	subss	%xmm6, %xmm4
	movss	36(%rdi), %xmm6
	mulss	%xmm0, %xmm6
	mulss	12(%rdi), %xmm0
	subss	%xmm6, %xmm4
	addss	%xmm5, %xmm0
	movss	4(%rdi), %xmm5
	mulss	%xmm4, %xmm5
	addss	%xmm5, %xmm0
	movss	%xmm0, 4(%rax,%r11,8)
	movq	48(%rdi), %rbx
	addq	$1, %r11
	movss	%xmm3, 48(%rdi)
	movq	%rbx, 56(%rdi)
	movss	%xmm4, 52(%rdi)
	cmpl	%r11d, %ecx
	jg	.L30
	movl	%r14d, %r11d
	xorl	%ebx, %ebx
	.p2align 4
	.p2align 3
.L31:
	leal	2182(%r11), %ebp
	leal	2525(%r11), %r12d
	andl	$32767, %ebp
	andl	$32767, %r12d
	movss	4(%rdx,%r12,8), %xmm0
	movss	(%rdx,%rbp,8), %xmm3
	movl	%r11d, %ebp
	subl	$1, %r11d
	andl	$32767, %ebp
	andl	$32767, %r11d
	unpcklps	%xmm0, %xmm3
	movq	(%rax,%rbx,8), %xmm0
	movaps	%xmm3, %xmm4
	mulps	%xmm1, %xmm4
	subps	%xmm4, %xmm0
	movq	%xmm0, %r12
	mulps	%xmm1, %xmm0
	addps	%xmm3, %xmm0
	movlps	%xmm0, (%rax,%rbx,8)
	addq	$1, %rbx
	movq	%r12, (%rdx,%rbp,8)
	cmpl	%ebx, %ecx
	jg	.L31
	movl	%r14d, %r11d
	xorl	%ebx, %ebx
	.p2align 4
	.p2align 3
.L32:
	leal	5215(%r11), %ebp
	leal	4722(%r11), %r12d
	andl	$32767, %ebp
	andl	$32767, %r12d
	movss	4(%rdx,%r12,8), %xmm0
	movss	(%rdx,%rbp,8), %xmm3
	leal	2525(%r11), %ebp
	subl	$1, %r11d
	andl	$32767, %ebp
	andl	$32767, %r11d
	unpcklps	%xmm0, %xmm3
	movq	(%rax,%rbx,8), %xmm0
	movaps	%xmm3, %xmm4
	mulps	%xmm1, %xmm4
	subps	%xmm4, %xmm0
	movq	%xmm0, %r12
	mulps	%xmm1, %xmm0
	addps	%xmm3, %xmm0
	movlps	%xmm0, (%rax,%rbx,8)
	addq	$1, %rbx
	movq	%r12, (%rdx,%rbp,8)
	cmpl	%ebx, %ecx
	jg	.L32
	movl	%r14d, %r11d
	xorl	%ebx, %ebx
	.p2align 4
	.p2align 3
.L33:
	leal	7888(%r11), %ebp
	leal	8062(%r11), %r12d
	andl	$32767, %ebp
	andl	$32767, %r12d
	movss	4(%rdx,%r12,8), %xmm0
	movss	(%rdx,%rbp,8), %xmm3
	leal	5215(%r11), %ebp
	subl	$1, %r11d
	andl	$32767, %ebp
	andl	$32767, %r11d
	unpcklps	%xmm0, %xmm3
	movq	(%rax,%rbx,8), %xmm0
	movaps	%xmm3, %xmm4
	mulps	%xmm1, %xmm4
	subps	%xmm4, %xmm0
	movq	%xmm0, %r12
	mulps	%xmm1, %xmm0
	addps	%xmm3, %xmm0
	movlps	%xmm0, (%rax,%rbx,8)
	addq	$1, %rbx
	movq	%r12, (%rdx,%rbp,8)
	cmpl	%ebx, %ecx
	jg	.L33
	movl	%r14d, %r11d
	xorl	%ebx, %ebx
	.p2align 4
	.p2align 3
.L34:
	leal	11492(%r11), %ebp
	leal	11272(%r11), %r12d
	andl	$32767, %ebp
	andl	$32767, %r12d
	movss	4(%rdx,%r12,8), %xmm0
	movss	(%rdx,%rbp,8), %xmm3
	leal	8062(%r11), %ebp
	subl	$1, %r11d
	andl	$32767, %ebp
	andl	$32767, %r11d
	unpcklps	%xmm0, %xmm3
	movq	(%rax,%rbx,8), %xmm0
	movaps	%xmm3, %xmm4
	mulps	%xmm1, %xmm4
	subps	%xmm4, %xmm0
	movq	%xmm0, %r12
	mulps	%xmm1, %xmm0
	addps	%xmm3, %xmm0
	movlps	%xmm0, (%rax,%rbx,8)
	addq	$1, %rbx
	movq	%r12, (%rdx,%rbp,8)
	cmpl	%ebx, %ecx
	jg	.L34
	movl	%r14d, %r11d
	xorl	%ebx, %ebx
	.p2align 6
	.p2align 4
	.p2align 3
.L35:
	movq	(%rax,%rbx,8), %r12
	leal	11492(%r11), %ebp
	addq	$1, %rbx
	subl	$1, %r11d
	andl	$32767, %ebp
	andl	$32767, %r11d
	movq	%r12, (%rdx,%rbp,8)
	cmpl	%ebx, %ecx
	jg	.L35
	cmpl	$1, -104(%rsp)
	jle	.L48
	movq	-24(%rsp), %r11
	leaq	4(%r13), %r12
	leaq	(%rax,%r11), %rbx
	movq	-112(%rsp), %r11
	addq	%r15, %r11
	cmpq	%r11, %rax
	setnb	%r11b
	cmpq	%rbx, %r15
	setnb	%bpl
	orl	%ebp, %r11d
	movq	%r15, %rbp
	subq	%r12, %rbp
	cmpq	$8, %rbp
	seta	%bpl
	testb	%bpl, %r11b
	je	.L48
	movq	-112(%rsp), %r11
	addq	%r13, %r11
	cmpq	%r11, %rax
	setnb	%r11b
	cmpq	%rbx, %r13
	setnb	%bl
	orb	%r11b, %bl
	je	.L48
	cmpl	$3, -104(%rsp)
	jle	.L49
	xorl	%r11d, %r11d
	.p2align 6
	.p2align 4
	.p2align 3
.L38:
	movups	(%rax,%r11,2), %xmm0
	movups	16(%rax,%r11,2), %xmm3
	movaps	%xmm0, %xmm4
	shufps	$221, %xmm3, %xmm0
	shufps	$136, %xmm3, %xmm4
	movups	%xmm4, 0(%r13,%r11)
	movups	%xmm0, (%r15,%r11)
	addq	$16, %r11
	cmpq	%r8, %r11
	jne	.L38
	movq	-80(%rsp), %r11
	movl	-56(%rsp), %ebx
	leaq	(%rax,%r11), %r12
	movq	-88(%rsp), %r11
	leaq	0(,%r11,4), %rbp
	leaq	(%r15,%rbp), %r11
	movq	%r11, -72(%rsp)
	leaq	0(%rbp,%r13), %r11
	movq	%r11, -64(%rsp)
	testb	$3, -100(%rsp)
	je	.L43
.L37:
	movl	-100(%rsp), %r11d
	subl	%ebx, %r11d
	cmpl	$1, %r11d
	je	.L40
	leaq	(%rax,%rbx,8), %rbp
	movq	0(%rbp), %xmm0
	movq	8(%rbp), %xmm3
	movaps	%xmm0, %xmm4
	unpcklps	%xmm3, %xmm0
	unpcklps	%xmm3, %xmm4
	shufps	$78, %xmm0, %xmm0
	movlps	%xmm4, 0(%r13,%rbx,4)
	movlps	%xmm0, (%r15,%rbx,4)
	testb	$1, %r11b
	je	.L43
	andl	$-2, %r11d
	leaq	(%r12,%r11,8), %r12
	salq	$2, %r11
	addq	%r11, -72(%rsp)
	addq	%r11, -64(%rsp)
.L40:
	movss	(%r12), %xmm0
	movq	-64(%rsp), %rbx
	movss	%xmm0, (%rbx)
	movq	-72(%rsp), %rbx
	movss	4(%r12), %xmm0
	movss	%xmm0, (%rbx)
.L43:
	movq	-96(%rsp), %rbx
	addq	%rbx, %r15
	addq	%rbx, %r13
.L18:
	movl	-116(%rsp), %r11d
	movq	-48(%rsp), %rbx
	subl	%ecx, %r14d
	andl	$32767, %r14d
	subl	%ecx, %r11d
	addq	%rbx, %rax
	andl	$2047, %r11d
	subl	%ecx, -52(%rsp)
	movl	%r11d, -116(%rsp)
	jne	.L44
	movl	%r14d, %r15d
.L16:
	movl	-116(%rsp), %eax
	movl	%r15d, 80(%rdi)
	movq	%rdx, 88(%rdi)
	movl	%eax, 64(%rdi)
	movq	-16(%rsp), %rax
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	movq	%rax, 72(%rdi)
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
.L48:
	.cfi_restore_state
	xorl	%r11d, %r11d
	.p2align 6
	.p2align 4
	.p2align 3
.L42:
	movss	(%rax,%r11,8), %xmm0
	movss	%xmm0, 0(%r13,%r11,4)
	movss	4(%rax,%r11,8), %xmm0
	movss	%xmm0, (%r15,%r11,4)
	addq	$1, %r11
	cmpl	%r11d, %ecx
	jg	.L42
	jmp	.L43
.L46:
	xorl	%r11d, %r11d
	.p2align 5
	.p2align 4
	.p2align 3
.L22:
	movss	(%r10,%r11,4), %xmm0
	addss	(%r9,%r11,4), %xmm0
	movss	%xmm0, (%rax,%r11,8)
	addq	$1, %r11
	cmpl	%r11d, %ecx
	jg	.L22
	jmp	.L23
.L49:
	movq	%r15, -72(%rsp)
	movq	%rax, %r12
	xorl	%ebx, %ebx
	movq	%r13, -64(%rsp)
	jmp	.L37
	.cfi_endproc
.LFE2309:
	.size	_ZN6Reverb7processEPPfS1_i, .-_ZN6Reverb7processEPPfS1_i
	.p2align 4
	.globl	_Z13SCReverb_nextP8SCReverbi
	.type	_Z13SCReverb_nextP8SCReverbi, @function
_Z13SCReverb_nextP8SCReverbi:
.LFB2339:
	.cfi_startproc
	subq	$56, %rsp
	.cfi_def_cfa_offset 64
	movq	72(%rdi), %rax
	addq	$104, %rdi
	movq	%fs:40, %rcx
	movq	%rcx, 40(%rsp)
	movl	%esi, %ecx
	leaq	16(%rsp), %rdx
	movq	%rsp, %rsi
	movdqu	(%rax), %xmm0
	movq	-24(%rdi), %rax
	movaps	%xmm0, (%rsp)
	movdqu	(%rax), %xmm0
	movaps	%xmm0, 16(%rsp)
	call	_ZN6Reverb7processEPPfS1_i
	movq	40(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L98
	addq	$56, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
.L98:
	.cfi_restore_state
	call	__stack_chk_fail@PLT
	.cfi_endproc
.LFE2339:
	.size	_Z13SCReverb_nextP8SCReverbi, .-_Z13SCReverb_nextP8SCReverbi
	.p2align 4
	.globl	api_version
	.type	api_version, @function
api_version:
.LFB2340:
	.cfi_startproc
	movl	$3, %eax
	ret
	.cfi_endproc
.LFE2340:
	.size	api_version, .-api_version
	.p2align 4
	.globl	server_type
	.type	server_type, @function
server_type:
.LFB2341:
	.cfi_startproc
	xorl	%eax, %eax
	ret
	.cfi_endproc
.LFE2341:
	.size	server_type, .-server_type
	.section	.rodata.str1.1
.LC14:
	.string	"SCReverb"
	.text
	.p2align 4
	.globl	load
	.type	load, @function
load:
.LFB2342:
	.cfi_startproc
	movq	%rdi, _ZL2ft(%rip)
	movq	48(%rdi), %rax
	xorl	%r8d, %r8d
	leaq	_Z13SCReverb_DtorP8SCReverb(%rip), %rcx
	leaq	_Z13SCReverb_CtorP8SCReverb(%rip), %rdx
	movl	$208, %esi
	leaq	.LC14(%rip), %rdi
	jmp	*%rax
	.cfi_endproc
.LFE2342:
	.size	load, .-load
	.local	_ZL2ft
	.comm	_ZL2ft,8,8
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC0:
	.long	1032463870
	.long	1032463870
	.long	1040852478
	.long	1040852478
	.align 16
.LC1:
	.long	1032463870
	.long	1032463870
	.long	1065353216
	.long	1065353216
	.align 16
.LC2:
	.long	-1080931025
	.long	-1080931025
	.long	1054038716
	.long	1054038716
	.set	.LC8,.LC0
	.set	.LC9,.LC1
	.set	.LC10,.LC2
	.section	.rodata.cst4,"aM",@progbits,4
	.align 4
.LC11:
	.long	1061158912
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC12:
	.long	1058139013
	.long	1058139013
	.align 8
.LC13:
	.long	-1086324736
	.long	-1086324736
	.ident	"GCC: (GNU) 14.2.1 20240910"
	.section	.note.GNU-stack,"",@progbits
